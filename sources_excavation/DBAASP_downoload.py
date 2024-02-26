from Bio import SeqIO
from tools.DBAASP_req import DBASSP_req

if __name__ == "__main__":
    sequences = list(SeqIO.parse('DBAASP.fasta', 'fasta'))
    with open("DBAASP.csv", 'w') as my_file:
        my_file.write(
            'ID\tName\tSequence\tseqLen\tsynthesisType\tcomplexity\ttargetGroups\ttargeObjects\tkingdom\tsubkingdom\tsource\tplasmid\tgene\n')
    for seq in sequences:
        try:
            req = DBASSP_req(id=seq.id, base_url='https://dbaasp.org/peptides/')
            rque_data = req.dowload_whole_data()
            str_row = ''
            for i in range(len(rque_data)):
                value = list(rque_data.values())[i]
                if i == len(rque_data) - 1:
                    str_row += f'{value}\n'
                else:
                    str_row += f'{value}\t'
            print(str_row)
            with open('DBAASP.csv', 'a') as my_file:
                my_file.write(str_row)
        except BaseException:
            with open('DBAASP.csv', 'a') as my_file:
                my_file.write(f'{seq.id}\t{seq.seq}\tJsonError!\t\n')