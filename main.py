from scripts.classify_ibd import current_ibd
from scripts.future_microbiome_state import print_bacteria, prepare_and_run_var, analyze_time_series_data, filtering

def current_ibd_analysis():
    current_ibd('data/LiuLloyd_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy_WNA.csv')
    current_ibd('data/LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy_WNA.csv.csv')
    current_ibd('data/Gevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized.csv')

def var_attempt():
    top = print_bacteria('data/bacteria.txt')
    filtering('data/Gevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized.csv', top)

    file_path = 'data/filtered.csv'
    datetime_column = 'Collection Week'
    prepared_data, var_results = prepare_and_run_var(file_path, datetime_column)

    if var_results:
        print(var_results.summary())
    analyze_time_series_data('data/filtered.csv')



if __name__ == "__main__":

    current_ibd_analysis()

    var_attempt()


