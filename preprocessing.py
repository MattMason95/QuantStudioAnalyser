## Import libraries
from typing import List, Optional, Union
from pathlib import Path
import pandas as pd
import numpy as np
from difflib import SequenceMatcher as SM
import os

class DataProcessor:
    '''
    A class handling the loading and processing of data from QuantStudio qPCR software.
    
    Attributes: 
    data: Input data to be processed. This can be provided directly as a pre-loaded dataframe or as a filepath location.
    
    '''

    def __init__(self, data: Union[pd.DataFrame,str,Path]):
        '''
        Initialise the DataProcessor with either a dataframe from memory or a filepath to a .csv file. 
        '''
        self.data = data  # Store the input data
        self.filepath = Path(data) if isinstance(data, (str, Path)) else None
        self.data = self._load_data()  # Process the data
        self._validate_data()

    def _load_data(self):
        '''
        Load user data from memory as a DataFrame or from a filepath. 
        '''
        if isinstance(self.data, pd.DataFrame):
            return self.data
        
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f'No file found at {self.filepath}.')
        
        file_extension = self.filepath.suffix.lower()
        if file_extension == '.csv':
            return pd.read_csv(self.filepath)
        elif file_extension in ['.xlsx', '.xls']:
            return pd.read_excel(self.filepath)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

    def _validate_data(self) -> None:
        '''
        Validate the format of the loaded data.
        '''
        # Basic DataFrame validation
        if not isinstance(self.data, pd.DataFrame):
            raise TypeError("Data not loaded as a DataFrame.")
        
        if self.data.empty:
            raise ValueError("Provided DataFrame is empty.")
            
        # File-specific validation
        if self.filepath is not None:
            # Check if file exists
            if not os.path.exists(self.filepath):
                raise FileNotFoundError(f"File no longer exists at path: {self.filepath}.")
            
            # Check if file is readable
            if not os.access(self.filepath, os.R_OK):
                raise PermissionError(f"No read permission for file: {self.filepath}.")
            
            # Check file size 
            file_size = os.path.getsize(self.filepath)
            limit = 500_000_000
            if file_size > limit:  # 500MB limit
                raise ValueError(f"File size ({file_size} bytes) exceeds limit ({limit/100_000}MB)")

    def input_function(self, availableGenes: str, availableConditions: List[str]) -> List:
        '''
        Wait for user to provide input values for the control gene(s) and the control condition(s).
        '''
        control_genes = str(input(f"Select control genes from {availableGenes}.")).split(',')
        
        control_conditions = str(input(f"Select control condition prefix from {availableConditions}."))
        
        invalid_condition = True
        while invalid_condition:
            # First check if any condition matches (case insensitive)
            matching_conditions = [condition for condition in availableConditions 
                                if condition.lower().startswith(control_conditions.lower())]
            
            if matching_conditions:
                # Extract the correct case prefix from the first matching condition
                # by taking the same number of characters as the input length
                correct_case_prefix = matching_conditions[0][:len(control_conditions)]
                control_conditions = correct_case_prefix
                print(f'Control condition identified: {control_conditions}')
                invalid_condition = False
            else:
                control_conditions = str(input(f"ERROR! Please submit a valid condition prefix from {availableConditions}."))
        
        return control_genes, control_conditions


    def parser(self) -> pd.DataFrame: 
        ## Parse data from input file
        ## Columns expected from the QuantStudio output

        expected_columns = ['Well Position','Sample Name','Target Name','CT','Ct Threshold','Tm1']

        ## Actual columns in the loaded data
        columns = self.data.columns.tolist()

        ## Are there any unexpected columns in the input data? 
        unexpected_columns = [x for x in columns if x not in expected_columns]

        print(f'Input data has columns: {columns}. There are {len(unexpected_columns)} unexpected column(s): {unexpected_columns}.')

        ## Declare qPCR targets found in the dataset
        targets = self.data['Target Name'].unique().tolist()
        print(f'Input data has {len(targets)} unique targets: {targets}.')

        ## Declare qPCR samples found in the dataset
        samples = self.data['Sample Name'].unique().tolist()
        print(f'Input data has {len(samples)} unique samples: {samples}.')
        
        ## Preserve the original data by making a copy for modification (memory constraints shouldn't be an issue for this type of data)
        parsed_data = self.data.copy()

        ## Drop any unexpected columns as these aren't
        if unexpected_columns: 
            parsed_data = parsed_data.drop(columns=unexpected_columns)

        ## Wait for user to provide input for the control genes and control columns (required for ddCT analysis)
        user_genes, user_conditions = self.input_function(availableGenes=targets, availableConditions=samples)

        ## To make user input easier and less error prone, incorporate some fuzzy string matching with sequence matcher
        ## Generate output container
        control_genes = []

        if len(user_genes) == 1:
            true_target = targets[np.argmax(SM(None,user_genes,x) for x in targets)]
            control_genes = true_target
      
        parsed_data['IsControlGene'] = parsed_data['Target Name'].apply(lambda x: 1 if x in control_genes else 0)
        parsed_data['IsControlCondition'] = parsed_data['Sample Name'].apply(lambda x: 1 if user_conditions in x else 0)

        return parsed_data

    
    def dataCleaning(self) -> pd.DataFrame:
        '''
        Clean data to detect and remove outliers
        '''
        ## Retrieve data from upstream parser function
        data = self.parser()
        
        ## Computing dCT for all conditions
        for name, grouped_data in data.groupby('IsControlCondition'):
            if name: 
                print('Control')
            else: 
                print('Target')
            
            mean_data = grouped_data.groupby(['Sample Name','Target Name']).mean().reset_index(drop=False)
            samples = mean_data['Sample Name'].unique()
            
            print(samples)
            return mean_data

            

        # control_data = data[data['IsControlCondition'] == 1]

        

        # for sample in samples:




        # return data

    # def reports(self,) -> dict:
    #     ## Generate reports from the outcomes of dataCleaning