## Import libraries
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from dataclasses import dataclass
from difflib import SequenceMatcher as SM
from typing import Dict, List, Optional, Union, Any

from file_logging import FileTracker

@dataclass
class FileProcessResult:
    ''' 
    Data class for main function. Returns an object containing the dataframe and additional metadata.
    '''
    dataframe: pd.DataFrame
    metadata: Dict[str, Any]

    def __getitem__(self, key):
        if key == 'dataframe':
            return self.dataframe
        elif key == 'metadata':
            return self.metadata
        else:
            raise KeyError(f'Invalid key: {key}.')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class DataProcessor:
    '''
    A class handling the loading and processing of data from QuantStudio qPCR software.
    
    Attributes: 
    data: Input data to be processed. This can be provided directly as a pre-loaded dataframe or as a filepath location.
    
    '''

    def __init__(self,
                 data: Union[pd.DataFrame,str,Path],
                 log_dir: str = 'logs'):
        '''
        Initialise the DataProcessor with either a dataframe from memory or a filepath to a .csv file. 
        '''
        self.file_tracker = FileTracker(log_dir=log_dir)
        self.filepath = Path(data) if isinstance(data, (str, Path)) else None

        if isinstance(data, pd.DataFrame):
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            file_info = {
                'operation':'dataframe_loaded',
                'file_info':None,
                'extra_details':None
            }
            self.data = FileProcessResult(dataframe=self.data, metadata=file_info)
        else: 
            self.data = self.LoadData()
        
        self._validate_data()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def LoadData(self):
        '''
        Load user data from memory as a DataFrame or from a filepath. 
        '''        
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f'No file found at {self.filepath}.')
        
        file_extension = self.filepath.suffix.lower()
        extra_info = {'file_extension':file_extension}

        if file_extension == '.csv':
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            ## Load data from .csv file 
            loaded_data = pd.read_csv(self.filepath)

            ## Generate unique file info, store metadata 
            file_info = self.file_tracker.track_operation(self.filepath, 'file_read', extra_info)
            
            ## Return dataframe and metadata in FileProcessResult object
            return FileProcessResult(dataframe=loaded_data,metadata=file_info)
        
        elif file_extension in ['.xlsx', '.xls']:
            print('\n##### Data loaded. #####\n'.center(shutil.get_terminal_size().columns))
            ## Load data direct from excel 
            loaded_data = pd.read_excel(self.filepath)

            ## Generate unique file information, store metadata
            file_info  = self.file_tracker.track_operation(self.filepath, 'file_read', extra_info)

            ## Return dataframe and metadata in FileProcessResult object
            return FileProcessResult(dataframe=loaded_data,metadata=file_info)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def _validate_data(self) -> None:
        '''
        Validate the format of the loaded data.
        '''
        # Basic DataFrame validation
        data = self.data.dataframe
        
        if not isinstance(data, pd.DataFrame):
            raise TypeError("Data not loaded as a DataFrame.")
        
        if data.empty:
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
        Optimise user input to find best match in the existing data. 
        Prompt re-entry if there are no good matches.

        Attributes:
        availableGenes - Targets extracted from the loaded data that could be selected as control genes.
        availableConditions - Sample identifiers extracted from the loaded data that could be selected as control conditions. 
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

    def Parser(self) -> FileProcessResult: 
        ## Parse data from input file
        ## Columns expected from the QuantStudio output

        print('\n##### Starting data parsing. #####\n'.center(shutil.get_terminal_size().columns))

        expected_columns = ['Well Position','Sample Name','Target Name','CT','Ct Threshold','Tm1']

        data = self.data

        # if data is not FileProcessResult:
        #     print('Data loaded incorrectly.')

        data = data.dataframe
        # metadata = data.metadata

        ## Actual columns in the loaded data
        columns = data.columns.tolist()

        ## Are there any unexpected columns in the input data? 
        unexpected_columns = [x for x in columns if x not in expected_columns]

        print(f'Input data has columns: {columns}. There are {len(unexpected_columns)} unexpected column(s): {unexpected_columns}.')

        ## Drop any unexpected columns as these aren't required downstream
        if unexpected_columns: 
            parsed_data = parsed_data.drop(columns=unexpected_columns)
            print(f'Removing unexpected columns: {unexpected_columns}')

        ## Declare qPCR targets found in the dataset
        targets = data['Target Name'].unique().tolist()
        print(f'Input data has {len(targets)} unique targets: {targets}.')

        ## Declare qPCR samples found in the dataset
        samples = data['Sample Name'].unique().tolist()
        print(f'Input data has {len(samples)} unique samples: {samples}.')
        
        ## Preserve the original data by making a copy for modification (memory constraints shouldn't be an issue for this type of data)
        parsed_data = data.copy()

        ## Wait for user to provide input for the control genes and control columns (required for ddCT analysis)
        user_genes, user_conditions = self.InputFunction(availableGenes=targets, availableConditions=samples)

        ## Designate control and target identities for genes and samples      
        parsed_data['TargetType'] = parsed_data['Target Name'].apply(lambda x: 'Control' if x in user_genes else 'Target')
        parsed_data['Condition'] = parsed_data['Sample Name'].apply(lambda x: 'Control' if user_conditions in x else 'Target')

        print('\n##### Completed data parsing. #####\n'.center(shutil.get_terminal_size().columns))
        
        # file_info = {
        #     'operation':'parsed',
        #     'N_targets':len(targets),
        #     'N_samples':len(samples)
        # }

        extra_info = {
            'N_targets':len(targets),
            'N_samples':len(samples)
        }
        
        file_info = self.file_tracker.track_operation(self.filepath, 'parsed', extra_info)

        return FileProcessResult(dataframe=parsed_data, metadata=file_info)

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
    def DataCleaning(self, data: pd.DataFrame):
        '''
        An embedded function to evaluate data for the detection and removal of outliers.
        '''
        ## Retrieve data from upstream function
        ## Get mean Tms for each target (condition independent)
        print('\n##### Starting data cleaning. #####\n'.center(shutil.get_terminal_size().columns))
    
        averages = {}
        variances = {}
        
        fig, axs = plt.subplots(4,1,figsize=(8,8))

        ax = axs.reshape(-1)
        
        for idx, (name, subset) in enumerate(data.groupby('Target Name')):
            i = idx + 1 
            averages[f'Tm1:{name}'] = np.median(subset['Tm1'])
            averages[f'CT:{name}'] = np.median(subset['CT'])
            variances[f'Tm1:{name}'] = subset['Tm1'].std()
            variances[f'CT:{name}'] = subset['CT'].std()
            
            ## Plotting of data
            ax[idx].boxplot(x=subset['Tm1'],vert=False,whis=0.75)
            ax[idx+i].boxplot(x=subset['CT'],vert=False,whis=0.75)
            ax.set_xlabel('Melting Temp.')
            ax.set_ylabel(name)

        ## Example indices for culling, following evaluate of medians/variances
        indices = [0,3,5,7] 

        print(f"Outlying data detected at indices: {indices}. Culling rows.")

        return indices, averages, variances

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def ddctCalculation(self, data=None, correction=False) -> FileProcessResult:
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
        
        if data is FileProcessResult:
            function_map = {
                'loaded':'DataProcessor()',
                'parsed':'Parser()'
            }
            if data.metadata['operation'] != 'parsed':
                raise ValueError(f'Data requires parsing. Please use Parser() function.')

        data = data.dataframe
        
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
            # indices, _,_ = self.DataCleaning(ddct_data)
            indices = [0,3,5,7]
            cleaned_data = ddct_data.drop(labels=indices)

            ## Get means of remaining technical replicates following outlier removal
            mean_ddCT = cleaned_data.groupby(['Sample Name','Target Name','Condition','TargetType']).mean().reset_index(drop=False)

            extra_info = {
                        'removed_data':indices
                        }
            
            file_info = self.file_tracker.track_operation(self.filepath, operation='ddct_calculation',extra_info=extra_info)

            return FileProcessResult(dataframe=mean_ddCT,metadata=file_info)

        ## Get means of all technical replicates
        mean_ddCT = ddct_data.groupby(['Sample Name','Target Name','Condition','TargetType']).mean().reset_index(drop=False)
        
        extra_info = {}

        file_info = self.file_tracker.track_operation(self.filepath, operation='ddct_calculation',extra_info=extra_info)

        ## Return the mean fold-change ddCT data for plotting
        return FileProcessResult(dataframe=mean_ddCT,metadata=file_info)