# importing necessary libraries 
import os
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split

# Path to the NEW364.csv dataset file
NEW364_PATH = r"C:\NEW364.csv"

# Training Parameters
VALIDATION_SPLIT = 0.2  # 20% of data for validation
BATCH_SIZE = 16
EPOCHS = 100
LEARNING_RATE = 1e-3
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Checkpoint and Logging
LOAD_CHECKPOINT_PATH = None

# Model Parameters
INPUT_DIM = 21      # 21 features per amino acid
HIDDEN_DIM = 128
OUTPUT_DIM = 8      # 8 classes for secondary structure
NUM_LAYERS = 2
# Max sequence length for padding.
MAX_SEQ_LEN = 700

# Data Loading
def load_csv_data(path, max_len=MAX_SEQ_LEN):
    """
    Loads and preprocesses data from a CSV file with 'input' and 'dssp8' columns.

    This function performs the following steps:
    1. Reads the CSV file using pandas.
    2. For each protein sequence:
       - One-hot encodes the amino acid sequence ('input') into a feature matrix.
       - One-hot encodes the 8-state DSSP labels ('dssp8') into a target matrix.
    3. Pads or truncates all sequences and labels to a fixed `max_len`.
    4. Creates a boolean mask to distinguish real residues from padding.
    """
    print(f"Loading and preprocessing dataset from: {path}")
    if not os.path.exists(path):
        print(f"Error: Dataset file not found at {path}")
        print("Please check the NEW364_PATH in the CONFIG section.")
        exit()

    # Define mappings for one-hot encoding
    # 20 standard amino acids + 'X' for unknown/other
    aa_map = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWYX')}
    # 8 DSSP classes. Note: The label for 'Coil' can vary (e.g., 'C', 'L', or ' ').
    # We assume 'C' for Coil here. Adjust if your dataset uses a different label.
    dssp8_labels = 'HBEGITS' + 'C' # H=Helix, B=Bridge, E=Strand, G=3-10 Helix, I=Pi-Helix, T=Turn, S=Bend, C=Coil
    dssp8_map = {label: i for i, label in enumerate(dssp8_labels)}

    df = pd.read_csv(path)
    num_proteins = len(df)

    # Initialize empty numpy arrays for the preprocessed data
    X = np.zeros((num_proteins, max_len, INPUT_DIM), dtype=np.float32)
    y = np.zeros((num_proteins, max_len, OUTPUT_DIM), dtype=np.float32)
    mask = np.zeros((num_proteins, max_len), dtype=np.bool_)

    for i, row in df.iterrows():
        seq = str(row['input'])
        dssp8 = str(row['dssp8'])

        # Ensure sequence and labels have the same length before processing
        min_len = min(len(seq), len(dssp8))
        # Determine the length to process, capped by max_len
        proc_len = min(min_len, max_len)

        # Set the mask for the actual length of the sequence
        mask[i, :proc_len] = True

        # Process input sequence (X) with one-hot encoding
        for j in range(proc_len):
            aa_idx = aa_map.get(seq[j].upper(), aa_map['X']) # Default to 'X' if amino acid is not standard
            X[i, j, aa_idx] = 1.0

        # Process target labels (y) with one-hot encoding
        for j in range(proc_len):
            label_idx = dssp8_map.get(dssp8[j].upper())
            if label_idx is not None:
                y[i, j, label_idx] = 1.0

    print("Dataset loaded and preprocessed successfully.")
    print(f"  - Input shape (X): {X.shape}")
    print(f"  - Target shape (y): {y.shape}")
    print(f"  - Mask shape (mask): {mask.shape}")
    return X, y, mask


# PyTorch Dataset

class ProteinDataset(Dataset):
    """PyTorch Dataset for protein sequence and structure data."""
    def __init__(self, X, y, mask):
        # Convert numpy arrays to PyTorch tensors
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32)
        self.mask = torch.tensor(mask, dtype=torch.bool)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.mask[idx]


# BiLSTM Model

class BiLSTM(nn.Module):
    """
    A Bidirectional LSTM model for sequence-to-sequence prediction.
    It processes a sequence of amino acid features and predicts a secondary
    structure class for each residue.
    """
    def __init__(self, input_dim=INPUT_DIM, hidden_dim=HIDDEN_DIM, output_dim=OUTPUT_DIM, num_layers=NUM_LAYERS):
        super().__init__()
        self.lstm = nn.LSTM(
            input_dim,
            hidden_dim,
            num_layers=num_layers,
            batch_first=True,       # Input format: (batch, seq_len, features)
            bidirectional=True      # Use information from past and future
        )
        # The output of the BiLSTM has 2 * hidden_dim features
        self.fc = nn.Linear(hidden_dim * 2, output_dim)

    def forward(self, x):
        # x shape: (batch, seq_len, input_dim)
        out, _ = self.lstm(x)      # out shape: (batch, seq_len, hidden_dim * 2)
        out = self.fc(out)         # out shape: (batch, seq_len, output_dim)
        return out


# Evaluation Function

def evaluate_model(model, dataloader, criterion, device):
    """Evaluates the model on a given dataset (e.g., validation set)."""
    model.eval()
    total_loss = 0
    total_correct = 0
    total_elements = 0

    with torch.no_grad():
        for batch_x, batch_y, batch_mask in dataloader:
            batch_x, batch_y, batch_mask = batch_x.to(device), batch_y.to(device), batch_mask.to(device)

            logits = model(batch_x)
            targets = batch_y.argmax(dim=-1)

            loss = criterion(logits.view(-1, OUTPUT_DIM), targets.view(-1))
            masked_loss = loss * batch_mask.view(-1).float()
            final_loss = masked_loss.sum() / batch_mask.sum()
            total_loss += final_loss.item()

            # Calculate accuracy (Q8)
            preds = logits.argmax(dim=-1)
            correct = (preds == targets) & batch_mask
            total_correct += correct.sum().item()
            total_elements += batch_mask.sum().item()

    avg_loss = total_loss / len(dataloader)
    accuracy = total_correct / total_elements if total_elements > 0 else 0
    return avg_loss, accuracy


# Training & Checkpointing

def save_checkpoint(state, filename="checkpoint.pth.tar"):
    """Saves model and training parameters at checkpoint."""
    print(" => Saving new best model...")
    torch.save(state, filename)

def train_and_evaluate(model, train_loader, val_loader, optimizer, criterion, epochs, device, save_dir, start_epoch=1, best_val_loss=float('inf')):
    """Handles the full training and validation loop, including checkpointing."""
    print(f"Starting training on {device} from epoch {start_epoch}...")

    for epoch in range(start_epoch, epochs + 1):
        # --- Training Phase ---
        model.train()
        total_train_loss = 0
        num_batches = len(train_loader)

        for i, (batch_x, batch_y, batch_mask) in enumerate(train_loader):
            batch_x, batch_y, batch_mask = batch_x.to(device), batch_y.to(device), batch_mask.to(device)
            optimizer.zero_grad()
            logits = model(batch_x)
            targets = batch_y.argmax(dim=-1)
            loss = criterion(logits.view(-1, OUTPUT_DIM), targets.view(-1))
            masked_loss = loss * batch_mask.view(-1).float()
            final_loss = masked_loss.sum() / batch_mask.sum()
            final_loss.backward()
            optimizer.step()

            total_train_loss += final_loss.item()

        avg_train_loss = total_train_loss / num_batches

        # Validation Phase
        val_loss, val_accuracy = evaluate_model(model, val_loader, criterion, device)

        print(
            f"Epoch {epoch}/{epochs} | "
            f"Train Loss: {avg_train_loss:.4f} | "
            f"Val Loss: {val_loss:.4f} | "
            f"Val Accuracy (Q8): {val_accuracy:.4f}"
        )

        # Checkpointing
        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss
            save_checkpoint({
                'epoch': epoch,
                'state_dict': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'best_val_loss': best_val_loss,
                'params': {
                    'batch_size': BATCH_SIZE,
                    'epochs': EPOCHS,
                    'lr': LEARNING_RATE,
                    'input_dim': INPUT_DIM,
                    'hidden_dim': HIDDEN_DIM,
                    'output_dim': OUTPUT_DIM,
                    'num_layers': NUM_LAYERS
                }
            }, filename=os.path.join(save_dir, 'best_model.pth.tar'))

    print("Training finished.")

# Inference
def load_model_for_inference(checkpoint_path, device="cpu"):
    """
    Loads a trained model from a checkpoint specifically for inference.
    It reconstructs the model architecture based on the saved parameters.
    """
    if not os.path.exists(checkpoint_path):
        print(f"Error: Checkpoint file not found at {checkpoint_path}")
        return None

    checkpoint = torch.load(checkpoint_path, map_location=device)
    
    # Recreate the model with the saved architecture parameters
    model_params = checkpoint['params']
    model = BiLSTM(
        input_dim=model_params['input_dim'],
        hidden_dim=model_params['hidden_dim'],
        output_dim=model_params['output_dim'],
        num_layers=model_params['num_layers']
    )
    
    # Load the trained weights
    model.load_state_dict(checkpoint['state_dict'])
    
    # Set the model to evaluation mode and move to the correct device
    model.to(device)
    model.eval()
    
    print("\nModel loaded successfully for inference.")
    return model

# Main Execution
if __name__ == "__main__":
    # Ensure the save directory exists for model checkpoints
    os.makedirs(SAVE_DIR, exist_ok=True)
    print("\n" + "="*60 + f"\n NEW RUN STARTED\n" + "="*60)
    
    # 1. Load data
    X, y, mask = load_csv_data(NEW364_PATH)

    # 2. Split data into training and validation sets
    indices = np.arange(X.shape[0])
    train_indices, val_indices = train_test_split(
        indices, test_size=VALIDATION_SPLIT, random_state=42, shuffle=True
    )
    
    train_dataset = ProteinDataset(X[train_indices], y[train_indices], mask[train_indices])
    val_dataset = ProteinDataset(X[val_indices], y[val_indices], mask[val_indices])

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)
    print(f"Data split: {len(train_dataset)} training samples, {len(val_dataset)} validation samples.")

    # 3. Initialize model, optimizer, and loss function
    model = BiLSTM().to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    criterion = nn.CrossEntropyLoss(reduction="none")
    print(f"Model initialized on {DEVICE}.")

    # 4. Load checkpoint if specified
    start_epoch = 1
    best_val_loss = float('inf')
    if LOAD_CHECKPOINT_PATH and os.path.exists(LOAD_CHECKPOINT_PATH):
        print(f"Loading checkpoint '{LOAD_CHECKPOINT_PATH}'")
        checkpoint = torch.load(LOAD_CHECKPOINT_PATH, map_location=DEVICE)
        model.load_state_dict(checkpoint['state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        start_epoch = checkpoint['epoch'] + 1
        best_val_loss = checkpoint['best_val_loss']
        print(f"Resuming training from epoch {start_epoch} with best validation loss {best_val_loss:.4f}")

    # 5. Start training and evaluation
    train_and_evaluate(
        model, train_loader, val_loader, optimizer, criterion,
        epochs=EPOCHS, device=DEVICE, save_dir=SAVE_DIR,
        start_epoch=start_epoch, best_val_loss=best_val_loss
    )

