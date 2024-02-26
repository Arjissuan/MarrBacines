import requests
import pandas as pd

class DBASSP_req:

    def __init__(self, id:str, seq:str = '', base_url:str ='https://dbaasp.org/peptides/'):
        self.id = id
        self.base_url = base_url
        self.seq = self.check_sequence(seq)

    def response(self) -> requests.models.Response:
        return requests.get(f'{self.base_url}{self.id}')

    def check_sequence(self, seq:str) -> str:
        if seq != '':
            return seq
        else:
            if self.response().status_code != 204:
                seqseq = self.response().json()['sequence']
                return seqseq
            else:
                return 'Error'

    def isin(self, value:str, vector:list[list[str]]) -> bool:
        for arr in vector:
            if value in arr:
                return True
        return False

    #nalezy dodac jeszcze targetGroups, targetObjects, i nwm potem sie spytam
    def data_collection(self):
        if self.response().status_code != 204:
            data = self.response().json()
            if data['synthesisType']['name'] == 'Ribosomal' or data['synthesisType']['name'] == 'NonRibosomal':
                if data['complexity']['name'] == 'Monomer':
                    func = lambda x: [x['source'], x['kingdom']['name'], x['subkingdom'], x['plasmid'], x['gene']]
                    sources_list = list(map(func, data['sourceGenes']))  # this one is working
                    return sources_list
                elif data['complexity']['name'] == 'Multimer' or data['complexity']['name'] == 'Multi-Peptide':
                    peptides = list(map(lambda x: x, data['monomers']))
                    #print(peptides[0]['sourceGenes'][0]['source'])
                    sourceGenes = list(map(lambda x: x['sourceGenes'], peptides))
                    sources_list = []
                    for monomer in sourceGenes:
                        for source in monomer:
                            warunek = self.isin(source['source'], sources_list)
                            #print(warunek)
                            if warunek is False:
                                sources_list.append([source['source'],
                                                     source['kingdom']['name'],
                                                     source['subkingdom'],
                                                     source['plasmid'],
                                                     source['gene']
                                                     ])
                    return sources_list

    def dowload_whole_data(self) -> dict:
        if self.response().status_code != 204:
            data = self.response().json()
            complxity = data['complexity']['name']
            synth_type = data['synthesisType']['name']
            db_dict = {'ID': self.id,
                       'name': data['name'],
                       'sequence': self.seq,
                       'seqLen': len(self.seq),
                       'synthesisType': synth_type,
                       'coplexity': complxity,
                       'targetGrups': list(map(lambda x: x['name'], data['targetGroups'])),
                       'targetObjects': list(map(lambda x: x['name'], data['targetObjects'])),
                       'kingdom': [],
                       'subkingdom': [],
                       'source': [],
                       'plasmid': [],
                       'gene': [],
                       }
            if synth_type != 'Synthetic':
                if complxity == 'Monomer':
                    sourceGenesfunc = lambda x: [x['source'],
                                                 x['kingdom']['name'],
                                                 x['subkingdom'],
                                                 x['plasmid'],
                                                 x['gene']
                                                 ]
                    sources_list = list(map(sourceGenesfunc, data['sourceGenes']))
                    db_dict['kingdom'] = list(map(lambda x: sources_list[x][1], range(len(sources_list))))
                    db_dict['subkingdom'] = list(map(lambda x: sources_list[x][2], range(len(sources_list))))
                    db_dict['source'] = list(map(lambda x: sources_list[x][0], range(len(sources_list))))
                    db_dict['plasmid'] = list(map(lambda x: sources_list[x][3], range(len(sources_list))))
                    db_dict['gene'] = list(map(lambda x: sources_list[x][4], range(len(sources_list))))
                else:
                    peptides = list(map(lambda x: x, data['monomers']))
                    # print(peptides[0]['sourceGenes'][0]['source'])
                    sourceGenes = list(map(lambda x: x['sourceGenes'][0], peptides))
                    # print(f'Source genes: {sourceGenes}')
                    for item in sourceGenes:
                        db_dict['kingdom'].append(item['kingdom']['name'])
                        db_dict['subkingdom'].append(item['subkingdom'])
                        db_dict['source'].append(item['source'])
                        db_dict['plasmid'].append(item['plasmid'])
                        db_dict['gene'].append(item['gene'])
            else:
                pass

            return db_dict


if __name__ == "__main__":
    test1 = DBASSP_req(id='DBAASPR_1')
    # print(test1.data_collection())
    print(test1.dowload_whole_data())
    test2 = DBASSP_req('DBAASPR_872')
    print(test2.dowload_whole_data())
    test3 = DBASSP_req("DBAASPR_21")
    # print(test3.data_collection())
    print(test3.dowload_whole_data())
    test4 = DBASSP_req('DBAASPS_6004')
    # print(test4.data_collection())
    print(test4.dowload_whole_data())
    test5 = DBASSP_req('DBAASPR_9561', 'QSHLSLCRYCCNCCRNKGCGYCCKF')
    # print(test5.data_collection())
    print(test5.dowload_whole_data())

    test6 = DBASSP_req(id='DBAASPN_21710')
    print(test6.dowload_whole_data())
