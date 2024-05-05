import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
from statsmodels.tsa.vector_ar.var_model import VAR

def analyze_time_series_data(filepath):
    df = pd.read_csv(filepath)
    df.dropna(inplace=True)
    subject_counts = df['Subject'].value_counts()
    repeating_subjects = subject_counts[subject_counts > 1]
    data_shape = df.shape

    print(f"DataFrame shape after preprocessing: {data_shape}")
    print("Counts of repeating 'Subject' entries:")
    print(repeating_subjects)

    analyze_subject_repetitions(filepath)

def analyze_subject_repetitions(filepath):
    df = pd.read_csv(filepath)
    df.dropna(inplace=True)

    print(f"DataFrame shape after preprocessing: {df.shape}")

    subject_counts = df['Subject'].value_counts()

    average_repeats = subject_counts.mean()

    print("Counts of repeating 'Subject' entries:")
    print(subject_counts)
    print(f"Average number of repeats per subject: {average_repeats:.2f}")

    return df, subject_counts, average_repeats


def print_bacteria(file_path):
    top_50ish = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            if len(parts) > 1:
                print(parts[0])
                top_50ish.append(parts[0])
            else:
                print("No colon found in:", line.strip())

    return top_50ish

def filter(file_path, bacteria_of_interest):
    df = pd.read_csv(file_path)
    filtered_df = df[bacteria_of_interest]
    filtered_df.to_csv('filtered.csv', index=False)

    print("Columns retained after filtering:", filtered_df.columns.tolist())

def prepare_and_run_var(filepath, datetime_col, numeric_cols_prefix='Bacteria', use_pca=True):
    df = pd.read_csv(filepath)
    df[datetime_col] = pd.to_datetime(df[datetime_col], errors='coerce')
    df.sort_values(by=datetime_col, inplace=True)
    df.set_index(datetime_col, inplace=True)

    if df.index.duplicated().any():
        df = df[~df.index.duplicated(keep='first')]

    if df.index.inferred_freq is None:
        df = df.asfreq('D')

    numeric_cols = [col for col in df.columns if numeric_cols_prefix in col]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

    df.dropna(subset=numeric_cols, inplace=True)
    if use_pca and df.shape[1] > 1 and df.shape[0] > 1:
        n_components = min(df.shape[0], df.shape[1], 5)
        pca = PCA(n_components=n_components)
        df[numeric_cols] = pca.fit_transform(df[numeric_cols])
        numeric_cols = ['PC' + str(i) for i in range(n_components)]
    elif use_pca:
        print("Not enough data for PCA. Proceeding without PCA.")

    df_numeric = df[numeric_cols]

    try:
        model = VAR(df_numeric)
        maxlags = min(len(df_numeric) // 3, 15)
        results = model.fit(maxlags=maxlags, ic='aic')
        print("VAR model fitted. Here's the summary:")
        print(results.summary())
    except Exception as e:
        print(f"An error occurred while fitting the VAR model: {e}")
        results = None

    return df_numeric, results