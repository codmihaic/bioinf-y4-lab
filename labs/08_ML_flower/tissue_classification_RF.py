#Use a GTEx Gene Expression Data Subset to implement and evaluate a Random Forest model to classify tissue types based on gene expression profiles
#The dataset should contain at least 100 genes randomly selected from GTEx data for tissue types for three distinct classes: Liver, Kidney and Brain
# Import necessary libraries
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
data = pd.read_csv("tissue_gene_expression.csv")  # Create a subset from the GTEx Gene expression dataset or use dataset from previous lab example
print("Dataset Preview:")
print(data.head())

# Split data into features (X) and target (y)
X = data.drop(columns=['Tissue_Type'])
y = data['Tissue_Type']

# Encode target labels if necessary (e.g., string to numeric)
from sklearn.preprocessing import LabelEncoder
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)
print("Classes:", label_encoder.classes_)

# Split into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y_encoded, test_size=0.2, random_state=42
)

# Initialize and train a Random Forest model
rf_model = RandomForestClassifier(random_state=42, n_estimators=100)
rf_model.fit(X_train, y_train)

# Predict on the test set
y_pred = rf_model.predict(X_test)

# Evaluate the model
print("Classification Report:")
print(classification_report(y_test, y_pred, target_names=label_encoder.classes_))

# Confusion Matrix
conf_matrix = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues',
            xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix for Tissue Type Prediction")
plt.show()

# Feature Importance
feature_importances = rf_model.feature_importances_
plt.figure(figsize=(10, 6))
plt.barh(X.columns, feature_importances)
plt.xlabel("Feature Importance")
plt.title("Random Forest Feature Importance")
plt.show()

#output classification report, confusion matrix and feature importance plot
#Use GridSearchCV to optimize parameters like n_estimators and max_depth
#Experiment with selecting the top 10 most important features and retraining the model
