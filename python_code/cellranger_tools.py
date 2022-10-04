import os


__author__ = 'quanch'
__discriptions__ = '''This script is a tool set that help processing the Single-Cell rna-seq SRA data from NCBI SRA database. '''
__title__ = 'Single-Cell rna-seq SRA data processing tool set'
__version__ = '0.1.0a'

class CellrangerTools:
    
    """_summary_
    """
    
    
    def __init__(self, SRRsLocation: str, workingDirectory: str, transcriptome: str, loadModules: bool = False):
        """_summary_

        Args:
            SRRsLocation (str): The location of the SraAccList.txt file, this file usually can get fron NCBI SRA database.
            workingDirectory (str): The directory that you want to download the SRA files and perform the analysis (like the SRA_home).
            transcriptome (str): the reference to the transcriptome file location
            loadModules (bool, optional): whether to load the module of sratoolkit. Defaults to False.
            downloadDirectory (str): The directory where the SRA files will be downloaded.
            allSras (Dict): The dictionary that contains all the SRA files id and the location of the SRA files.
            SRRs (List): The list of the SRA files id.
        """

        self.SRRsLocation = SRRsLocation
        self.SRRs = self.__get_srrs_list__()
        self.workingDirectory = workingDirectory
        self.downloadDirectory = workingDirectory + '/download'
        self.transcriptome = transcriptome
        self.allSras = None
        if loadModules:
            os.system('module load Bioinformatics; module load sratoolkit')
        
   
    def __get_srrs_list__(self):
        """This function is used to get the SRRs list from the SraAccList.txt file.

        Returns:
            List: The list of the SRA files id
        """
        with open(self.SRRsLocation, 'r') as f:
            srrs = f.read().splitlines()
        return srrs
    
    def prefetch(self):
        """This function is used to download the SRA files from NCBI SRA database using the sratoolkit.
        it will download the SRA files to the download directory and create a _fastq directory for the fastq-dump.

        Returns:
            self: return the CellrangerTools object for more operations.
        """
        os.system('cd ' + self.downloadDirectory)
        os.system('prefetch --option-file ' + self.SRRsLocation+' -O'+self.downloadDirectory)
        allSras = {}
        for srr in os.listdir(self.downloadDirectory):
            if srr in self.SRRs:
                sra_location = self.downloadDirectory + '/' + srr 
                os.makedirs(self.downloadDirectory + '/' + srr + '/' + '_fastq')
                allSras[srr] = sra_location
        self.allSras = allSras
        return self

    def fastq_dump(self):
        """This function is used to dump the single-cell rna-seq fastq files from the SRA files using the sratoolkit.
        it will dump the fastq files to the _fastq directory and change their name to the cellranger format.

        Returns:
            self: return the CellrangerTools object for more operations.
        """
        if self.allSras is None:
            self.__get_allSras__()

        for srr,path in self.allSras.items():
            sra_path = path + '/' + srr + '.sra'
            print(f'Performing fastq-dump on {srr}')
            os.system('fastq-dump --split-files ' + sra_path + ' --outdir ' + path + '/_fastq')
            self.__rename_fastq__(directory = path + '/' + '_fastq')
            print(f'Finish fastq-dump on {srr}.')
        return self
    
    def __get_allSras__(self):
        """This function is used to get the allSras dictionary from the download directory if exists.
        
        Raises:
            Exception: if the download directory is not exists, it will raise an exception.
            run the prefetch function first to generate the allSras dictionary.or the download directory.

        Returns:
            None: this function will return nothing, but it will set the self.allSras attribute.
        """
        if self.allSras is not None:
            return None
        downloadDir = os.listdir(self.downloadDirectory)
        if len(downloadDir) > 0:
            allSras = {}
            for srr in downloadDir:
                if srr in self.SRRs:
                    sra_location = self.downloadDirectory + '/' + srr
                    allSras[srr] = sra_location
                else:
                    raise Exception('Please run prefetch first.')
            self.allSras = allSras
        else:
            raise Exception('Please run prefetch first.')
        
    def __rename_fastq__(self, directory: str = None):
        map_dic = {'1':'I1', '2':'R1', '3':'R2'}
        for fastq in os.listdir(directory):
            SRR, index = fastq.split('.')[0].split('_')
            newFileName = f'{SRR}_S1_L001_{map_dic[index]}_001.fastq'
            old_path = directory+'/'+fastq
            new_path = directory+'/'+newFileName
            os.rename(old_path,new_path)
            print(f'{fastq} is renamed to {newFileName}')

    def cellranger_count(self,ids:list=None, transcriptome:str=None, cores:int=None, mem:int=None):
        """this function is used to perform the cellranger count function on the fastq files.

        Args:
            ids (List, optional): the id list get form SraAccList. Defaults to None. if None, it will use the self.SRRs attribute.
            transcriptome (str, optional): the reference to the transcriptome file location. Defaults to None. if None, it will use the self.transcriptome attribute.
            cores (int, optional): number of cores to use for cellranger. Defaults to None.
            mem (int, optional): number of memory to use for cellranger. Defaults to None.
        """
        
        if transcriptome is None:
            transcriptome = self.transcriptome
        if ids is None:
            ids = self.SRRs
        self.__get_allSras__()

        for id in ids:
            cwd = self.allSras[id]
            os.chdir(cwd)

            fastqs = cwd+"/_fastq"
            command = f'cellranger count --id={id} \
                        --transcriptome={transcriptome} \
                        --fastqs={fastqs} \
                        --sample={id}'
            if cores:
                command += f' --localcores={cores}'
            if mem:
                command += f' --localmem={mem}'

            os.system(command)
            os.rename(f'{cwd}/{id}', f'{cwd}/CountResult')

        
    @staticmethod
    def remove_all_files(directory: str):
        os.chdir(directory)
        os.system('rm -rf *')




# def get_parser():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-n', '--NameList', type=str, help='the SraAccList.txt file location')


# def main():
#     parser = get_parser()



if __name__=='__main__':
    pass
    #main()