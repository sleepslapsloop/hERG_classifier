import joblib
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    roc_auc_score,
)

X_train = pd.read_pickle("../Data/TestTrain/X_train.pkl")
y_train = pd.read_pickle("../Data/TestTrain/y_train.pkl")
X_test = pd.read_pickle("../Data/TestTrain/X_test.pkl")
y_test = pd.read_pickle("../Data/TestTrain/y_test.pkl")

# 1. Initialize the Random Forest
# n_jobs=-1 uses all available CPU cores to speed up training
rf_model = RandomForestClassifier(
    n_estimators=500,
    max_depth=30,  # Restricts trees from growing infinitely deep
    min_samples_leaf=2,  # Forces leaves to contain at least 2 molecules (stops memorization)
    max_features="sqrt",  # Limits the number of bits each tree can look at once
    class_weight="balanced",
    random_state=42,
    n_jobs=-1,
)

# 2. Train the model
print("Training the Random Forest model (this might take a minute)...")
rf_model.fit(X_train, y_train)

# 3. Predict on the test set
y_pred = rf_model.predict(X_test)
y_proba = rf_model.predict_proba(X_test)[:, 1]  # Probabilities for the active class

# 4. Print performance metrics
print("\n================ MODEL PERFORMANCE ================")
print(f"Accuracy:  {accuracy_score(y_test, y_pred):.4f}")
print(f"ROC-AUC:   {roc_auc_score(y_test, y_proba):.4f}")
print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# Define the filename you want to save it as
model_filename = "../Data/Models/herg_classifier_01.joblib"

# Save the trained model to disk
joblib.dump(rf_model, model_filename)

print(f"Model successfully saved to {model_filename}")

# Predict on the training set
y_train_pred = rf_model.predict(X_train)
y_train_proba = rf_model.predict_proba(X_train)[:, 1]

print("--- Overfitting Check ---")
print(
    f"Train Accuracy: {accuracy_score(y_train, y_train_pred):.4f} | Test Accuracy: {accuracy_score(y_test, y_pred):.4f}"
)
print(
    f"Train ROC-AUC:  {roc_auc_score(y_train, y_train_proba):.4f} | Test ROC-AUC:  {roc_auc_score(y_test, y_proba):.4f}"
)
