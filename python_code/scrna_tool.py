import os


__author__ = 'quanch'
__discriptions__ = '''This script is a tool set that help processing the Single-Cell rna-seq SRA data from NCBI SRA database. '''
__title__ = 'Single-Cell rna-seq SRA data processing tool set'
__version__ = '0.1.0a'

class ScRNAAnalysis:
    def __init__(self,SraAccListPath, SRAHomePath,transcriptomePath):
        self.__SraAccListPath = SraAccListPath
        self.__SRRList = None
        self.__SraHomePath = SRAHomePath
        self.__transcriptomePath = transcriptomePath
        self.__SRRLocationDict = {}
    
    #region property
    @property
    def SraAccListPath(self):
        return self.__SraAccListPath
    @property
    def SRAHomePath(self):
        return self.__SraHomePath
    @property
    def transcriptomePath(self):
        return self.__transcriptomePath
    @property
    def SRRList(self):
        return self.__SRRList
    @property
    def SRRLocationDict(self):
        return self.__SRRLocationDict
    #endregion

    #region methods
    def _get_SRRs_list(self):
        """This function is used to get the SRRs list from the SraAccList.txt file.

        Returns:
            List: The list of the SRA files id
        """
        SRRList = []
        with open(self.SraAccListPath, 'r') as f:
            for line in f:
                if line.startswith('SRR'):
                    SRRList.append(line.strip())
        self.__SRRList = SRRList
        
        for SRR in SRRList:
            self.__SRRLocationDict[SRR] = {}
        return self
    
    def _check_local_file(self):
        """This funciton is used to check if the SRR file is already downloaded
        in the SRAHomePath.

        Returns:
        List of the SRRs that need to be prefetch
        """
        if self.__SRRList is None:
            self._get_SRRs_list()
        SRRsNeedToPrefetch = self.__SRRList.copy()
        SRRList = self.__SRRList.copy()
        SRAHomePath = self.__SraHomePath

        for SRR in SRRList:
            SRR_folder_path = os.path.join(SRAHomePath, SRR)
            if os.path.exists(SRR_folder_path):
                # check if SRR.sra exists
                SRR_sra_path = os.path.join(SRR_folder_path, SRR + '.sra')
                if os.path.exists(SRR_sra_path):
                    SRRsNeedToPrefetch.remove(SRR)
                    self.__SRRLocationDict[SRR]['SRR_dir'] = SRR_folder_path
                else:
                    print('SRR.sra file does not exist in {}, prefect later.'.format(SRR_folder_path))
        return SRRsNeedToPrefetch
    
    def prefetch_SRRs(self):
        """This function is used to prefetch the SRRs that need to be downloaded.
        """
        SRRsNeedToPrefetch = self._check_local_file()

        for SRR in SRRsNeedToPrefetch:
            print('prefetching {}...'.format(SRR))
            os.system('prefetch {}'.format(SRR))
            print('prefetching {} finished.'.format(SRR))
            # update the SRRLocationDict
            self.__SRRLocationDict[SRR]['SRR_dir'] = os.path.join(self.__SraHomePath, SRR)
        return self
    
    def fastq_dump(self):
        """This function is used to dump the single-cell rna-seq fastq files from the SRA files using the sratoolkit.
        it will dump the fastq files to the _fastq directory and change their name to the cellranger format.

        Returns:
            self: return the CellrangerTools object for more operations.
        """
        if len(self.__SRRLocationDict) == 0:
            self.prefetch_SRRs()
        for SRRid, SRRpath in self.__SRRLocationDict.items():
            SRRFilePath = os.path.join(SRRpath, SRRid + '.sra')
            fastqDir = os.path.join(SRRpath, '/_fastq')
            print('dumping {}...'.format(SRRid))
            os.system(f'fastq-dump {SRRFilePath}\
                --split-files --gzip --outdir {fastqDir}')
            self.__rename_fastq__(fastqDir)
            # update the SRRLocationDict
            self.__SRRLocationDict[SRRid]['fastq_dir'] = fastqDir
            print('dumping {} finished.'.format(SRRid))

        return self

    def parallel_fastq_dump(self,threads=1):
        os.chdir(self.__SraHomePath)
        if len(self.__SRRLocationDict) == 0:
            self.prefetch_SRRs()
        for SRRid, SRRDict in self.__SRRLocationDict.items():
            SRRpath = SRRDict['SRR_dir']
            os.chdir(SRRpath)
            FastqDir = SRRpath+'/fastqs'
            print('dumping {}...'.format(SRRid))
            os.system(f'parallel-fastq-dump --sra-id {SRRid}.sra --threads {threads} --outdir fastqs/ --split-files --gzip')
            self.__rename_fastq__(FastqDir)
            # update the SRRLocationDict
            self.__SRRLocationDict[SRRid]['fastq_dir'] = FastqDir
            print('dumping {} finished.'.format(SRRid))

        return self

    def __rename_fastq__(self, directory: str = None):
        map_dic = {'1':'I1', '2':'R1', '3':'R2'}
        for fastq in os.listdir(directory):
            SRR, index = fastq.split('.')[0].split('_')
            newFileName = f'{SRR}_S1_L001_{map_dic[index]}_001.fastq.gzip'
            old_path = directory+'/'+fastq
            new_path = directory+'/'+newFileName
            os.rename(old_path,new_path)
            print(f'{fastq} is renamed to {newFileName}')
        return self

    def cellranger_count(self,ids:list=None, transcriptome:str=None, cores:int=None, mem:int=None):
        """this function is used to perform the cellranger count function on the fastq files.

        Args:
            ids (List, optional): the id list get form SraAccList. Defaults to None. if None, it will use the self.SRRs attribute.
            transcriptome (str, optional): the reference to the transcriptome file location. Defaults to None. if None, it will use the self.transcriptome attribute.
            cores (int, optional): number of cores to use for cellranger. Defaults to None.
            mem (int, optional): number of memory to use for cellranger. Defaults to None.
        """
        
        if transcriptome is None:
            self.__transcriptomePath
        if ids is None:
            ids = self.__SRRList

        for id in ids:
            cwd = self.__SRRLocationDict[id]['SRR_dir']
            os.chdir(cwd)
            cellranger_path = '/media/evan/Linux_data/Bioinformatics/cellranger/cellranger-7.0.1/bin/cellranger'

            fastqs = cwd+"/fastqs"
            command = f'{cellranger_path} count --id={id} \
                        --transcriptome={transcriptome} \
                        --fastqs={fastqs} \
                        --sample={id}'
            if cores:
                command += f' --localcores={cores}'
            if mem:
                command += f' --localmem={mem}'

            os.system(command)
            #os.system(cellranger_path)
            os.rename(f'{cwd}/{id}', f'{cwd}/cellranger_{id}')
            self.__SRRLocationDict[id]['cellranger_dir'] = f'{cwd}/cellranger_{id}'
    
    def build_SRRdict(self):
        """This function is used to build the SRRLocationDict attribute.
        """
        
        # build the SRR_dir and fastq_dir
        for SRR in self.__SRRList:
            SRR_folder_path = os.path.join(self.__SraHomePath, SRR)
            if os.path.exists(SRR_folder_path):
                self.__SRRLocationDict[SRR]['SRR_dir'] = SRR_folder_path
                self.__SRRLocationDict[SRR]['fastq_dir'] = os.path.join(SRR_folder_path, '_fastq')
            else:
                print('SRR.sra file does not exist in {}, prefect later.'.format(SRR_folder_path))
        # build the cellranger_dir
        for SRR, SRRDict in self.__SRRLocationDict.items():
            SRR_folder_path = SRRDict['SRR_dir']
            cellranger_path = os.path.join(SRR_folder_path, 'cellranger_'+SRR)
            if os.path.exists(cellranger_path):
                self.__SRRLocationDict[SRR]['cellranger_dir'] = cellranger_path
            else:
                print('cellranger file does not exist in {}, cellranger later.'.format(SRR_folder_path))
        return self
    

