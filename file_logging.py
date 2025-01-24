## Import dependencies
import os
import json
import logging
import hashlib
import pandas as pd
from pathlib import Path
from typing import Union, Dict, Any, Optional
# from preprocessing import FileProcessResult


class FileTracker:
    '''
    Class for handling data and event logging from the main DataProcessor function.
    '''    
    def __init__(self,
                 log_dir: str = 'logs',
                 log_level: int = logging.INFO):     
        ## Make a log directory - if exists, leave unaltered
        os.makedirs(log_dir, exist_ok=True)
        
        ## Set logger and logging level
        self.logger = logging.getLogger('FileTracker')
        self.logger.setLevel(log_level)
        
        ## Set default file tracking log file at the log_dir path
        file_handler = logging.FileHandler(os.path.join(log_dir, 'file_tracking.log'))
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
        self.logger.addHandler(file_handler)
        
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        self.logger.addHandler(console_handler)
    
    def register_file(self,
                      filepath: Union[Path, str],
                      ) -> Dict[str, Any]:

        try:
            with open(filepath, 'rb') as f:
                file_hash = hashlib.md5(f.read()).hexdigest()
            
            return {
                'filename': os.path.basename(filepath),
                'size_bytes': os.path.getsize(filepath),
                'full_path': os.path.abspath(filepath),
                'hash': file_hash,
                'modified_time': os.path.getmtime(filepath)
            }
        except Exception as e:
            self.logger.error(f"Error generating file signature: {e}")
            return {}
        
            
    def track_operation(self, 
                        filepath: Union[Path, str],
                        operation: str, 
                        extra_info: Optional[Dict[str, Any]] = None):
               
        try:
            file_signature = self.register_file(filepath)
            
            log_entry = {
                'operation': operation,
                'file_info': file_signature,
                'extra_details': extra_info or {}
            }
            
            self.logger.info(json.dumps(log_entry, indent=2))

            meta_data = {
                'operation': operation,
                'file_info': file_signature,
            }

            return meta_data
        
        except Exception as e:
            self.logger.error(f"Tracking error for {filepath}: {e}")