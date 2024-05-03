import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

def prepare_data(df):
    # Convert timestamp into a datetime object and sort data
    df['Collection Week'] = pd.to_datetime(df['Collection Week'])
    df.sort_values(by=['Subject', 'Collection Week'], inplace=True)

    # Create lag features or differences that might indicate changes over time
    df['Previous_Diagnosis'] = df.groupby('Subject')['Diagnosis'].shift(1)
    df.fillna({'Previous_Diagnosis': 'No Diagnosis'}, inplace=True)  # Fill NA for first entry of each patient

    return df

def future_ibd():
    # Load data
    df = pd.read_csv('LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy.csv')

    # Assuming the first 5 columns are 'Sample', 'Diagnosis', 'Sex', 'Age', 'Sample Type', 'Timestamp'
    feature_columns = df.columns[5:-1].copy()  # Ensure to make a copy of the columns
    X = df[feature_columns].copy()  # Ensure to make a copy of the DataFrame slice
    y = df['Diagnosis'].copy()  # Ensure to make a copy of the Series

    # One-hot encode categorical variables
    X_encoded = pd.get_dummies(X, columns=[col for col in X.columns if X[col].dtype == 'object'])

    # Handle missing values
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Check and handle infinite values
    X_encoded.replace([np.inf, -np.inf], np.nan, inplace=True)
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_encoded)

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    # Train classifier
    classifier = RandomForestClassifier(random_state=42)
    classifier.fit(X_train, y_train)

    # Evaluate classifier
    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print("###########FOR FUTURE PREDICTIONS#############")
    print(f"Accuracy: {accuracy:.2f}")

    # Cross-validation for more robust evaluation
    scores = cross_val_score(classifier, X_scaled, y, cv=5)
    print(f"Average Cross-Validated Accuracy: {scores.mean():.2f}")

    # Get feature importances
    feature_importances = classifier.feature_importances_
    importance_df = pd.DataFrame(feature_importances, index=X_encoded.columns, columns=['Importance'])
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    print("Top 20 features by importance:")
    for index, row in importance_df.head(21).iterrows():
        print(f"{index}: {row['Importance']:.6f}")

def current_ibd():
    # Load data
    df = pd.read_csv('LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy.csv', low_memory=False)

    # Select columns for microbiome data
    microbiome_data = df.columns[5:55]
    X = df[microbiome_data]
    y = df['Diagnosis']

    # Identify and one-hot encode categorical columns
    X_encoded = pd.get_dummies(X, columns=[col for col in X.columns if X[col].dtype == 'object'])

    # Handle missing values
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Check and handle infinite values
    X_encoded.replace([np.inf, -np.inf], np.nan, inplace=True)
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_encoded)

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    # Train classifier
    classifier = RandomForestClassifier(random_state=42)
    classifier.fit(X_train, y_train)

    # Evaluate classifier
    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print("###########FOR CURRENT PREDICTIONS#############")
    print(f"Accuracy: {accuracy:.2f}")

    # Cross-validation
    scores = cross_val_score(classifier, X_scaled, y, cv=5)
    print(f"Average Cross-Validated Accuracy: {scores.mean():.2f}")

    # Get feature importances
    feature_importances = classifier.feature_importances_
    importance_df = pd.DataFrame(feature_importances, index=X_encoded.columns, columns=['Importance'])
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    print("Top 20 features by importance:")
    for index, row in importance_df.head(21).iterrows():
        print(f"{index}: {row['Importance']:.6f}")

if __name__ == '__main__':
    #current_ibd()
    future_ibd()