## Import libraries
from typing import List, Optional, Union
from pathlib import Path
import pandas as pd
import numpy as np
from difflib import SequenceMatcher as SM
import os
import shutil
import matplotlib.pyplot as plt

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
        self.data = self.LoadData()  # Process the data
        self._validate_data()
        self.manifest = False

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def LoadData(self):
        '''
        Load user data from memory as a DataFrame or from a filepath. 
        '''
        if isinstance(self.data, pd.DataFrame):
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            return self.data
        
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f'No file found at {self.filepath}.')
        
        file_extension = self.filepath.suffix.lower()
        if file_extension == '.csv':
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            return pd.read_csv(self.filepath)
        elif file_extension in ['.xlsx', '.xls']:
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            return pd.read_excel(self.filepath)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def InputFunction(self, availableGenes: str, availableConditions: List[str]) -> List:
        '''
        Wait for user to provide input values for the control gene(s) and the control condition(s).
        '''
        print(f'\n##### User Inputs. #####\n'.center(shutil.get_terminal_size().columns))

        control_genes = str(input(f"Select control genes from {availableGenes}.")).split(',')
        
        if len(control_genes) > 1:
                control_genes = str(input(f"ERROR! Too many genes submitted. Please select one control gene from {availableGenes}."))
        
        invalid_gene = True
        while invalid_gene: 
                      
            control_genes = ''.join(control_genes )

            ratios = [SM(None,''.join(control_genes),x).ratio() for x in availableGenes]

            ## If the maximum ratio is greater than 0.5 (i.e. a match better than 50%)
            if np.max(ratios) > .33:
                print(f'Control gene(s): {control_genes}.')
                ## Get index of highest match and corresponding gene from availableGenes
                target = availableGenes[np.argmax(ratios)]
                ## Replace the user-provided target with the closest matching target
                control_genes = target

                invalid_gene = False
            else: 
                control_genes = str(input(f"ERROR! No match found. Please submit a valid control gene from {availableGenes}."))
        
        control_conditions = str(input(f"Select control condition prefix from {availableConditions}."))

        invalid_condition = True
        while invalid_condition:
            # First check if any condition matches (case insensitive)
            matching_conditions = [condition for condition in availableConditions 
                                if condition.lower().startswith(control_conditions.lower())]
            
            if matching_conditions:
                print(f'Control condition(s): {control_conditions}.\n')
                # Extract the correct case prefix from the first matching condition
                # by taking the same number of characters as the input length
                correct_case_prefix = matching_conditions[0][:len(control_conditions)]
                control_conditions = correct_case_prefix
                
                invalid_condition = False
            else:
                control_conditions = str(input(f"ERROR! Please submit a valid condition prefix from {availableConditions}."))

        print('.'*shutil.get_terminal_size().columns,'\n')
        print(f'Control gene identified: {control_genes}.')
        print(f'Control condition identified: {control_conditions}.')
        return control_genes, control_conditions

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def Parser(self) -> pd.DataFrame: 
        ## Parse data from input file
        ## Columns expected from the QuantStudio output

        print('\n##### Starting data parsing. #####\n'.center(shutil.get_terminal_size().columns))

        expected_columns = ['Well Position','Sample Name','Target Name','CT','Ct Threshold','Tm1']

        ## Actual columns in the loaded data
        columns = self.data.columns.tolist()

        ## Are there any unexpected columns in the input data? 
        unexpected_columns = [x for x in columns if x not in expected_columns]

        print(f'Input data has columns: {columns}. There are {len(unexpected_columns)} unexpected column(s): {unexpected_columns}.')

        ## Drop any unexpected columns as these aren't required downstream
        if unexpected_columns: 
            parsed_data = parsed_data.drop(columns=unexpected_columns)
            print(f'Removing unexpected columns: {unexpected_columns}')

        ## Declare qPCR targets found in the dataset
        targets = self.data['Target Name'].unique().tolist()
        print(f'Input data has {len(targets)} unique targets: {targets}.')

        ## Declare qPCR samples found in the dataset
        samples = self.data['Sample Name'].unique().tolist()
        print(f'Input data has {len(samples)} unique samples: {samples}.')
        
        ## Preserve the original data by making a copy for modification (memory constraints shouldn't be an issue for this type of data)
        parsed_data = self.data.copy()

        ## Wait for user to provide input for the control genes and control columns (required for ddCT analysis)
        user_genes, user_conditions = self.InputFunction(availableGenes=targets, availableConditions=samples)

        ## Designate control and target identities for genes and samples      
        parsed_data['TargetType'] = parsed_data['Target Name'].apply(lambda x: 'Control' if x in user_genes else 'Target')
        parsed_data['Condition'] = parsed_data['Sample Name'].apply(lambda x: 'Control' if user_conditions in x else 'Target')

        print('\n##### Completed data parsing. #####\n'.center(shutil.get_terminal_size().columns))
        
        self.manifest = True

        return parsed_data

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    def reports(func) -> None:
        ## Decorator for generating reports from the DataCleaning function
        def wrapper(self, data):
            ## Run DataCleaning within wrapper
            outputs  = func(self, data)

            print(f'Compiling reports from {func.__name__}.')
            
            indices, averages, variances = outputs

            outlier_info = data.loc[indices]

            outlier_info['Tm1DeviationFromMedian'] = outlier_info.apply(lambda row: abs(row['Tm1'] - averages[f'Tm1:{row["Target Name"]}'])/variances[f'Tm1:{row["Target Name"]}'], axis=1)
            outlier_info['CtDeviationFromMedian'] = outlier_info.apply(lambda row: abs(row['CT'] - averages[f'CT:{row["Target Name"]}'])/variances[f'CT:{row["Target Name"]}'], axis=1)

            print('Stored outlier reports in DataProcessor().reports_data.')
            self.reports_data = outlier_info
            
            return outputs

        return wrapper

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    @reports
    def DataCleaning(self, data: pd.DataFrame) -> pd.DataFrame:
        '''
        An embedded function to evaluate data for the detection and removal of outliers.
        '''
        ## Retrieve data from upstream function
        ## Get mean Tms for each target (condition independent)
        print('\n##### Starting data cleaning. #####\n'.center(shutil.get_terminal_size().columns))
    
        averages = {}
        variances = {}
        
        fig, axs = plt.subplots(4,1,figsize=(8,8))

        for ax, (name, subset) in zip(axs.reshape(-1), data.groupby('Target Name')):
            averages[f'Tm1:{name}'] = np.median(subset['Tm1'])
            averages[f'CT:{name}'] = np.median(subset['CT'])
            variances[f'Tm1:{name}'] = subset['Tm1'].std()
            variances[f'CT:{name}'] = subset['CT'].std()
            
            ## Plotting of data
            ax.boxplot(x=subset['Tm1'],vert=False,whis=0.75)
            ax.boxplot(x=subset['CT'],vert=False,whis=0.75)
            ax.set_xlabel('Melting Temp.')
            ax.set_ylabel(name)

        ## Example indices for culling, following evaluate of medians/variances
        indices = [0,3,5,7] 

        print(f"Outlying data detected at indices: {indices}. Culling rows.")

        return indices, averages, variances

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def ddctCalculation(self, data=None, correction=False) -> pd.DataFrame:
        '''
        Clean data to detect and remove outliers
        
        Attributes: 
        - Correction [False / True] (default = False). Set true to remove outliers from data.  
        '''
        if data is None:         
            ## Retrieve data from upstream parser function
            if not hasattr(self, 'Parser') or not callable(getattr(self, 'Parser', None)):
                raise AttributeError('No Parser() method available to fetch data.')

            data = self.Parser()
        
        if not self.manifest:
            raise ValueError('Data has not been parsed. Please use the Parser() function first.')


        ## Prepare output dataframe for dCT data
        dct_data = pd.DataFrame()

        ## Computing dCT for all conditions
        for _, grouped_conditions in data.groupby('Condition'):            
            ## Could calculate means of technical replicates at this stage
            mean_data = grouped_conditions#.groupby(['Sample Name','Target Name','Condition','TargetType']).mean().reset_index(drop=False)
            
            ## Get mean CT  value for the housekeeping gene(s)
            hk_gene_mean = mean_data[mean_data['TargetType'] == 'Control']['CT'].mean()

            ## Calculate dCT by subtracting CT values from mean housekeeping gene(s) CT value
            mean_data['dCT'] = mean_data['CT'].apply(lambda x: hk_gene_mean - x)
            
            ## Concatenate output container
            dct_data = pd.concat([dct_data,mean_data])

        ## Prepare output dataframe for ddCT data
        ddct_data = pd.DataFrame()

        ## Computing ddCT for all targets
        for _, grouped_targets in dct_data.groupby('Target Name'):
            ## Get mean dCT for the control condition
            target_mean = grouped_targets[grouped_targets['Condition'] == 'Control']['dCT'].mean()
            
            ## Calculate ddCT by subtracting dCT values from the mean control condition dCT
            grouped_targets['ddCT'] = grouped_targets['dCT'].apply(lambda x: target_mean - x)

            ## Concatenate output container
            ddct_data = pd.concat([ddct_data,grouped_targets])

        ## Calculate fold-change with the 2^(-ddCT) method
        ddct_data['2^(-ddCT)'] = ddct_data['ddCT'].apply(lambda x: 2**(-x))

        if correction: 
            '''
            Pass data to DataCleaning() function to assess distributions of CTs/TMs 
            DataCleaning returns list of indices that meet criteria for culling
            Drop these indices from ddct_data and calculate means from technical replicates
            '''
            indices, _,_ = self.DataCleaning(ddct_data)
            cleaned_data = ddct_data.drop(labels=indices)

            ## Get means of remaining technical replicates following outlier removal
            mean_ddCT = cleaned_data.groupby(['Sample Name','Target Name','Condition','TargetType']).mean().reset_index(drop=False)
            return mean_ddCT

        ## Get means of all technical replicates
        mean_ddCT = ddct_data.groupby(['Sample Name','Target Name','Condition','TargetType']).mean().reset_index(drop=False)
        
        ## Return the mean fold-change ddCT data for plotting
        return mean_ddCT