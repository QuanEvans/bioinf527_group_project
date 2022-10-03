import os


__author__ = 'quanch'
__discriptions__ = '''This script is a tool set that help processing the Single-Cell rna-seq SRA data from NCBI SRA database. '''
__title__ = 'Single-Cell rna-seq SRA data processing tool set'
__version__ = '0.1.0a'

class CellrangerTools:
    
    """_summary_
    """
    
    
    def __init__(self, SRRsLocation: str, workingDirectory: str, cellrangerPath: str):
        """_summary_

        Args:
            SRRs (str): the SraAccList.txt file location
        """
        self.SRRs = self.get_srrs_list()
        self.workingDirectory = workingDirectory
        self.downloadDirectory = workingDirectory + '/download'
        self.cellrangerPath = cellrangerPath
        self.allSras = None
    
        
    def get_srrs_list(self):
        with open(self.SRRsLocation, 'r') as f:
            srrs = f.read().splitlines()
        return srrs
    
    def prefect_srrs(self):
        os.system('cd ' + self.downloadDirectory)
        os.system('prefetch --option-file ' + self.SRRsLocation)
        allSras = []
        for srr in os.listdir(self.downloadDirectory):
            if srr in self.SRRs:
                sra_location = self.downloadDirectory + '/' + srr + '/' + srr + '.sra'
                allSras.append(sra_location)
        self.allSras = allSras
        return self
        
    
    def fastq_dump(self):
        os.system('cd ' + self.workingDirectory)
        os.system('fastq-dump --split-file *.sra')
        return self
        
    
    def __rename_fastq__(self, directory: str = None):
        os.system('cd ' + directory)
        os.system('rename \'s/_1.fastq/_I1_001.fastq/\' *.fastq')
        os.system('rename \'s/_2.fastq.gz/_R1_001.fastq.gz/\' *.fastq.gz')
        os.system('rename \'s/_3.fastq.gz/_R2_001.fastq.gz/\' *.fastq.gz')




# def get_parser():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-n', '--NameList', type=str, help='the SraAccList.txt file location')


# def main():
#     parser = get_parser()



if __name__=='__main__':
    pass
    #main()