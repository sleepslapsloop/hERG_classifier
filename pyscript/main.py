import joblib
import numpy as np
from ProcessMol import compute_descriptors, embed_3d, smiles_to_mol

# 1. Load your trained Random Forest model
model_path = "../Data/Models/herg_classifier_01.joblib"
model = joblib.load(model_path)

# 2. Input SMILES string (Example: a known blocker)
input_smiles = "CN(C)c1ccc(/N=N/c2ccccn2)c(O)c1"

print("Step 1: Parsing SMILES and generating 3D conformer...")

# 3. Run the processing pipeline
mol_2d = smiles_to_mol(input_smiles)
mol_3d = embed_3d(mol_2d)
features = compute_descriptors(mol_3d)

# 4. Check if the molecule was successfully processed
if features is None:
    print(
        "❌ Error: Pipeline failed. Verify SMILES validity or 3D conformer generation."
    )
else:
    print("Step 2: Running model inference...")

    # Scikit-learn models expect a 2D array (shape: [1, n_features]) for single-sample predictions
    features_reshaped = features.reshape(1, -1)

    # 5. Get the binary prediction and class probabilities
    prediction = model.predict(features_reshaped)[0]
    probabilities = model.predict_proba(features_reshaped)[0]

    prob_safe = probabilities[0]
    prob_inhibitor = probabilities[1]

    # 6. Print the results clearly
    print("\n================= hERG PREDICTION RESULTS =================")
    print(f"SMILES:          {input_smiles}")
    print(
        f"Verdict:         {r'🚨 hERG INHIBITOR (Cardiotoxic Risk)' if prediction == 1 else r'✅ NON-INHIBITOR (Safe)'}"
    )
    print("-----------------------------------------------------------")
    print(f"Probability of being Safe (Class 0):      {prob_safe:.2%}")
    print(f"Probability of being Inhibitor (Class 1): {prob_inhibitor:.2%}")
    print("===========================================================")
