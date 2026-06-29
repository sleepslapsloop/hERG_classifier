import pandas as pd
from sklearn.model_selection import train_test_split

df = pd.read_pickle("../Data/Mol3D/hERG_ChEMBL_Mol3D_descriptors.pkl")

X = df.drop(columns=["Smiles", "pIC50", "Mol3D"])
y = (df["pIC50"] >= 5.5).astype(int)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

X_train.to_pickle("../Data/TestTrain/X_train.pkl")
X_test.to_pickle("../Data/TestTrain/X_test.pkl")
y_train.to_pickle("../Data/TestTrain/y_train.pkl")
y_test.to_pickle("../Data/TestTrain/y_test.pkl")

print(
    f"X_train: {X_train.shape}\nX_test: {X_test.shape}\ny_train: {y_train.shape}\ny_test: {y_test.shape}"
)
