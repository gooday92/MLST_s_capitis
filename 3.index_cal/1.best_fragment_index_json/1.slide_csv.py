from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import json
from icecream import ic

# 计算每个基因每个可变框的分辨能力，输出每个克隆和总体的index，每个克隆的主要等位基因，该处的等位基因编号不等于最终编号



def align_df(path):
    alignment = AlignIO.read(file_path, "fasta")
    genes_list = []
    for record in alignment:
        genes_list.append(record.id)
    df = pd.DataFrame(alignment, index=genes_list)
    return df


def save_json(result, path):
    json_str = json.dumps(result, indent=4)
    with open(path, "w") as json_file:
        json_file.write(json_str)

def best_slide(df_alin, ds_group, fra_length, gene):
    slide_list = []
    gene_length = len(df_alin.columns)
    result_index = dict()
    # 找到存在差异的突变位点，slide的起始位点从这个列表里遍历
    for col in range(0, gene_length-fra_length+1):
        column = df_alin.iloc[:,col].value_counts()
        if len(column) >= 2:
            slide_list.append(col)

    # slide滑动框计算index，和其他字典
    for col_start in slide_list:
        df_slide = df_alin.iloc[:, col_start:col_start+fra_length]
        sum_ = df_slide.isin(["-"]).sum().sum()
        if sum_:
            pass
        else:
            # print(df_slide)
            result_index[col_start] = dict()
            type_dict = dict()
            for strain, row in df_slide.iterrows():
                type_seq = "".join(row)
                strain = strain.split(";")[0]
                type_dict[strain] = type_seq
            ds = pd.Series(type_dict)
            df_groups = pd.DataFrame()
            df_groups["type"] = ds
            df_groups["group"] = ds_group
            # 编号不同的等位基因
            l = list(df_groups.type.value_counts().index)
            allels_type = dict()
            for i, seq in enumerate(l):
                allels_type[seq] = i + 1
            # print(allels_type)
            # 构建结果df
            ds_all = df_groups.type.value_counts().values
            result_index[col_start]["All"] = inx_count(ds_all)
            for cha in ["B", "C", "D", "E", "F", "L", "A", "Basal"]:
                df_g = df_groups[df_groups.group == cha]
                ds2 = df_g.type
                ds_value = ds2.value_counts()
                # print(gene, ds_value)
                index = inx_count(ds_value.values)
                result_index[col_start][cha] = index
                # 结果中增加每个克隆最多的等位基因编号和占比
                allel_id, percent = cal_percent(ds_value, allels_type)
                gene_percent = f"{allel_id}({percent})"
                label = f"{cha}_num"
                result_index[col_start][label] = gene_percent
    if len(result_index):
        df_result = pd.DataFrame(result_index).T
    else:
        df_result = pd.DataFrame()
    return df_result


def cal_percent(ds_value, allels_type):
    allel_id, percent = 0, 0
    if len(ds_value) == 0:
        allel_id, percent = 0, 0
    else:
        allel_id = allels_type[ds_value.index[0]]
        percent = round(ds_value[0]/ sum(ds_value), 4)
    return allel_id, percent

def inx_count(ds):
    sum = ds.sum()
    sum_numerator = 0
    for count in ds:
        sum_numerator += count*(count -1)
    index = 1 - (sum_numerator / (sum *(sum -1)))
    return index



path_input = r"D:\python\3.MLST_SC\5.MLST\1.pangenome_analysis\2_genes_align"
path_output = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\1.best_fragment_index_json\slid_json"
genes_list = ['abcA', 'accA', 'accD_2~~~accD~~~accD_1', 'ackA_2~~~ackA_1~~~ackA', 'acoA', 'acsA2', 'acsA_1~~~acsA_2~~~acsA', 'acuA', 'acuC', 'ada', 'addA~~~addA_1', 'addB', 'adk', 'aes_2~~~aes_1', 'aes_2~~~aes_1~~~aes_3', 'agrA', 'ahpC', 'alaS', 'aldA', 'ald_2~~~ald~~~ald_1', 'alr2', 'alsS', 'amoA_2~~~amoA_1~~~amoA', 'ampA', 'ampC_2~~~ampC~~~ampC_1', 'ampS_1~~~ampS_2~~~ampS', 'ansA_2~~~ansA~~~ansA_1', 'apaH~~~apaH_2', 'apbA', 'apt', 'araB~~~araB_1~~~araB_2', 'arcA_1~~~arcA_2~~~arcA', 'arcB_3~~~arcB_2~~~arcB~~~arcB_1', 'arcC_1~~~arcC~~~arcC_2', 'arcD_1~~~arcD_2~~~arcD_3~~~arcD', 'argE_1', 'argE_2', 'argG_1~~~argG_2~~~argG', 'argH_2~~~argH_1~~~argH', 'argS_1~~~argS_2~~~argS', 'arg~~~arg_2~~~arg_1', 'arlR', 'arlS_2~~~arlS_1~~~arlS', 'arnC_1', 'aroA', 'aroB', 'aroC', 'asd', 'asnC', 'asp1_2~~~asp1~~~asp1_1', 'asp2~~~asp2_1~~~asp2_2', 'aspS', 'atl_1~~~atl~~~atl_2', 'atl~~~atl_2', 'atpA_2', 'atpA_2~~~atpA_1', 'atpB_1', 'atpB_2', 'atpD', 'atpE', 'atpG', 'aur_1~~~aur_2~~~aur', 'ausA_2~~~ausA_1~~~ausA', 'azlC~~~azlC_1~~~azlC_2', 'azo1_1~~~azo1_2~~~azo1', 'azoR_1~~~azoR~~~azoR_2', 'bcp', 'bcr~~~bcr_1', 'bdbD_2~~~bdbD_1', 'betA', 'betT', 'bfmB', 'bfmBAB', 'bglG_2~~~bglG~~~bglG_1', 'bicA~~~dauA', 'birA', 'bshC', 'budA1', 'butA_1~~~butA~~~butA_2', 'caiA_1~~~caiA_2~~~caiA', 'camS', 'capA', 'capB', 'capC', 'capD', 'capD_2~~~ggt', 'carB', 'cbiO_1', 'cbiO_1~~~cbiO_2', 'ccpA', 'cdaR', 'cdsA_2~~~cdsA_1~~~cdsA', 'cidC', 'cinA_2~~~cinA_1~~~cinA', 'citB', 'citC_2~~~citC', 'citM_2~~~citM~~~citM_1', 'citZ', 'clpB', 'clpP', 'clpQ', 'clpX', 'cls1_1~~~cls1~~~cls1_2', 'cls2_2~~~cls2_1~~~cls2', 'cmpR~~~cysL', 'cntI', 'coaBC', 'cobW', 'codY', 'coiA', 'comEA_1~~~comEA_2~~~comEA', 'comEC_1~~~comEC_2~~~comEC', 'comFC~~~comFC_1~~~comFC_2', 'comK~~~comK_2~~~comK_1', 'copA_3~~~copA_2~~~copA~~~copA_1', 'corA_1~~~corA_2~~~corA', 'cpdA', 'crr', 'crtI_1~~~crtI_2~~~crtI', 'crtN_1~~~crtN_2~~~crtN_3~~~crtN', 'csbC~~~csbC_2~~~csbC_1', 'cstR', 'ctaA', 'ctaB', 'ctrA', 'ctsR', 'cudT_1~~~cudT_2~~~cudT', 'cvfB_1~~~cvfB_2~~~cvfB', 'cydB~~~cydB_1~~~cydB_2', 'cysC', 'cysE', 'cysG', 'cysH~~~cysH_2~~~cysH_1', 'cysI', 'cysJ~~~cysJ_2~~~cysJ_1', 'cysK', 'cysM', 'cysS', 'czrB', 'dapA', 'dapB_1~~~dapB_2~~~dapB', 'dapD', 'dcuA', 'ddl', 'def', 'def2', 'deoC2_1~~~deoC2_2~~~deoC2', 'deoD', 'deoR_2', 'dhoM_2~~~dhoM_1~~~dhoM', 'dinG', 'dinP~~~dinB', 'disA', 'div1b~~~divIB', 'divIC', 'dltA~~~dltA_1~~~dltA_2', 'dltB', 'dltD', 'dnaA', 'dnaB_2~~~dnaB~~~dnaB_1', 'dnaD', 'dnaE_2~~~dnaE', 'dnaI', 'dnaN', 'dnaX', 'dps', 'drm', 'dtd', 'ebpS', 'emrA', 'engA', 'engB', 'eno', 'est', 'eutD', 'exuT', 'ezrA', 'fabB', 'fabD', 'fabG_2', 'fabG_2~~~fabG_1~~~fabG_3', 'fabH_1', 'fabI', 'fabZ', 'fapR', 'fbaA', 'fbpA', 'fdaB', 'fdhD_1~~~fdhD_2~~~fdhD', 'fecB_2', 'fecB_2~~~fecB_1', 'fecD', 'fhs_1~~~fhs_2~~~fhs', 'fhuA_1~~~fpuD~~~fhuA_2~~~cbiO_1~~~btuD', 'fhuA_2~~~fhuA_1', 'fhuB~~~fhuB_1~~~fhuB_2', 'fhuD', 'fmhA_2~~~fmhA_1~~~fmhA', 'fmhB', 'fmt', 'fmtC~~~fmtC_2~~~fmtC_1', 'fni', 'folC', 'folD', 'folE', 'folK', 'folP', 'frr', 'ftn_1~~~ftn_2~~~ftn', 'ftsA', 'ftsK~~~ftsK_1~~~ftsK_2', 'ftsL', 'ftsW', 'ftsY', 'ftsZ', 'fur', 'galM_1~~~galM_2~~~galM', 'galU', 'gap', 'gapB~~~gapB_1~~~gapB_2', 'gatA', 'gatB', 'gbsA_1~~~gbsA_2', 'gbsA_2~~~aldA_1~~~gbsA_1', 'gcaD', 'gerCC_1~~~gerCC_2~~~gerCC', 'gidA', 'gidB', 'gid~~~gid_1~~~gid_2', 'glcT~~~glcT_2~~~glcT_1', 'glcU', 'glmM', 'glmS', 'glnA~~~glnA_2~~~glnA_1', 'glnQ_1~~~glnQ~~~glnQ_2', 'glpD', 'glpK', 'glpP', 'glpQ', 'gltB2_2~~~gltB2_1~~~gltB2', 'gltB_1~~~gltB~~~gltB_2', 'gltC_2~~~gltC~~~gltC_1', 'gltD_2~~~gltD_1~~~gltD', 'gltS_1~~~gltS_2~~~gltS', 'gltT_2~~~gltT~~~gltT_1', 'gltX', 'glxK_1~~~glxK~~~glxK_2', 'glyA', 'gmk', 'gnd', 'gntK_1', 'gntP', 'gntR_1~~~gntR_2~~~gntR', 'gph', 'gppA', 'gpsA_1~~~gpsA_2~~~gpsA', 'gpxA2', 'graR', 'graS~~~graS_1~~~graS_2', 'greA', 'groEL', 'gsaB', 'gtfA_1', 'gtfA~~~pglH', 'guaA', 'guaB', 'gudB_1', 'gudB_2~~~gudB_1~~~ldh_1', 'gyrA', 'gyrB', 'hemA', 'hemC_2~~~hemC~~~hemC_1', 'hflX~~~hflX_1~~~hflX_2', 'hisS', 'hmrA~~~hmrA_2~~~hmrA_1', 'holB', 'hprK', 'hslU', 'hssR', 'hssS_2~~~hssS~~~hssS_1', 'htsC', 'hutG', 'hutH_1~~~hutH_2~~~hutH', 'hutI', 'hutU_3~~~hutU_2~~~hutU_1~~~hutU', 'hxlA_1', 'icaC_1~~~icaC_2~~~icaC_3~~~icaC', 'ileS_2~~~ileS_1~~~ileS', 'ilvA_1', 'ilvB_1~~~ilvB~~~ilvB_2', 'ilvC_1~~~ilvC~~~ilvC_2', 'ilvD_1~~~ilvD_2~~~ilvD', 'ilvE', 'infB', 'infC', 'ipdC', 'isaB_1~~~isaB_2', 'isaB_2', 'iscU', 'isdB~~~harA', 'isdE', 'isdF~~~btuC~~~isdF_2~~~isdF_1', 'ispE', 'kbl', 'kdpE~~~kdpD~~~rcsC', 'kdtB', 'kgd', 'ksgA_1~~~ksgA_2~~~ksgA', 'lacA', 'lacB', 'lacC_2~~~lacC_1~~~lacC', 'lacD', 'lacG_1~~~lacG_2~~~lacG', 'larB', 'larC_1~~~larC_2~~~larC', 'larE', 'ldcC', 'ldh2', 'ldh~~~ldh_2~~~ldh_1', 'leuA_1~~~leuA_2~~~leuA_3~~~leuA', 'leuB_2~~~leuB_1~~~leuB', 'leuC_1~~~leuC~~~leuC_2', 'leuD_1~~~leuD~~~leuD_2', 'leuS', 'lexA_2~~~lexA_1~~~lexA', 'lgt', 'lig', 'lipA', 'lip_2~~~lip_1', 'lip_2~~~lip_1~~~lip_3~~~lip', 'lpdA', 'lplA2', 'lrgB_1~~~lrgB_2', 'lrgB_2', 'lspA~~~lspA_1', 'luxS', 'lysA', 'lysC', 'lysP', 'lytA', 'lytH', 'maeB', 'malA_3~~~malA_1~~~malA_2~~~malA', 'manA1_1~~~manA1~~~manA1_2', 'mazG_2~~~mazG_3', 'mcsB', 'mdrP', 'mdtD~~~mdtL', 'mecA2', 'menH', 'metB_1~~~metB_2~~~metB', 'metE', 'metH_2~~~metH_1~~~metH', 'metK', 'metS_2~~~metS_1~~~metS', 'metXA~~~metX_1~~~metX_2~~~metX', 'mfd', 'mgtE', 'mhpD_2~~~mhpD~~~mhpD_1', 'miaA_2~~~miaA', 'miaB_2~~~miaB_1~~~miaB', 'mnaA', 'mnhA~~~mnhA_1~~~mnhA_2', 'mnhB_1', 'mnhB_1~~~mnhB_2', 'mnhD_1', 'mnhD_1~~~mnhD_2', 'mnhE_1', 'mnhE_1~~~mnhE_2', 'mnhG_2~~~mnhG_1', 'mnh_A1', 'mntA', 'mntB', 'mntC', 'moaA', 'moaB', 'moaC', 'moaE', 'mobA_2~~~mobA_1~~~mobA', 'mobB', 'modA_2~~~modA_1~~~modA', 'modB_1~~~modB_2~~~modB', 'moeA_2~~~moeA_1~~~moeA', 'mqo1', 'mqo2_1', 'mraW', 'mraY_1~~~mraY_2~~~mraY', 'mraZ', 'mreC', 'msbA_2~~~msbA~~~msbA_1', 'msmX_2~~~btuD_2~~~btuD~~~potA_1~~~potA_2', 'msmX_2~~~msmX~~~msmX_1', 'msrB', 'msrR~~~msrR_2~~~msrR_1', 'mta~~~tipA', 'mtlA_1', 'mtlA_3~~~mtlA_2', 'mtlD_2~~~mtlD_1~~~mtlD', 'murA', 'murB~~~murB_2~~~murB_1', 'murC', 'murD', 'murE', 'murE2', 'murF', 'murG', 'murI', 'murP', 'murQ_1~~~murQ_2~~~murQ', 'murZ', 'mutL~~~mutL_1~~~mutL_2', 'mutS', 'mutY', 'mvaA', 'mvaD', 'mvaK1_1~~~mvaK1_2~~~mvaK1', 'mvaK2_3~~~mvaK2_2~~~mvaK2_1~~~mvaK2', 'mvaS', 'nagA', 'nagB_2~~~nagB_1~~~nagB', 'narG_2~~~narG_1~~~narG', 'narI', 'narJ_1~~~narJ_2~~~narJ', 'narY_1~~~narY_2~~~narY', 'nasF_1~~~nasF_2', 'nasF_3~~~nasF_2~~~nasF_1', 'ndhF', 'ndk', 'nhaP2', 'nifS', 'nifU_1~~~nifU_2~~~nifU', 'nikD_1~~~nikD~~~nikD_2', 'nikE', 'nirB~~~nirB_2~~~nirB_1', 'noc', 'norA_1~~~norA_2~~~norA', 'npr', 'nrdD_2~~~nrdD_1~~~nrdD', 'nrdE', 'nrdF', 'nrdG', 'nrdR', 'nreB_2~~~nreB_1~~~nreB', 'nrgA_2~~~nrgA_1~~~nrgA', 'nrmA', 'nth', 'nudF', 'nupC', 'nupC2_2~~~nupC2_1~~~nupC2', 'nusA', 'nusG', 'obgE', 'odhB', 'ohr', 'opp2B_2~~~opp2B_1~~~opp2B', 'opp2C_1~~~opp2C_2~~~opp2C', 'opuCA_2~~~opuCA_1', 'opuCB_1', 'opuCC_2~~~opuCC_1', 'opuCD_1~~~opuCD_2', 'osmC_2~~~osmC_1', 'oxaI_1', 'oxaI_2', 'pabA', 'pabB', 'pabC', 'panB', 'panC_1~~~panC_2~~~panC', 'panE', 'panK', 'pap', 'papS', 'parB', 'parC', 'parE', 'pbp2', 'pbp4~~~pbp4_1~~~pbp4_2', 'pbpA', 'pbuX_1~~~pbuX~~~pbuX_2', 'pckA~~~pckA_2~~~pckA_1', 'pcp_1~~~pcp_2~~~pcp', 'pcrA', 'pcrB~~~pcrB_2~~~pcrB_1', 'pdhA', 'pdhB', 'pdhC', 'pdhD', 'pdp_2~~~pdp_1~~~pdp', 'pdxS', 'pepA1', 'pepA2_1~~~pepA2_2~~~pepA2', 'pepT_1~~~pepT_2~~~pepT', 'pepX_1~~~pepX_2~~~pepX', 'pfkA', 'pfkB', 'pfpI~~~pfpI_1', 'pgi', 'pgk', 'pgm', 'pheA', 'pheS', 'pheT', 'pheT2_1~~~pheT2~~~pheT2_2', 'phoA', 'phoU', 'phrB_1~~~phrB_2~~~phrB', 'pitA', 'pknB', 'plsC', 'plsX', 'plsY', 'pmdD', 'pmi', 'pnbA_2~~~pnbA_1~~~pnbA', 'pnpA', 'polA', 'polC', 'porA', 'porB', 'potA~~~potA_2~~~potA_1', 'potB', 'potC_1~~~potC_2~~~potC', 'potD_2~~~potD_1~~~potD', 'ppdK_1~~~ppdK~~~ppdK_2', 'ppi', 'ppk', 'ppnK', 'prfA', 'prfC_1~~~prfC_2~~~prfC', 'priA', 'proC', 'proP_2~~~proP_1~~~proP', 'proS', 'prs', 'pstA', 'pstB', 'pstC', 'pstS~~~pstS_2', 'ptaA_1~~~ptaA~~~ptaA_2', 'pth', 'ptsA_1~~~ptsA_2~~~ptsA', 'purB', 'purC', 'purE', 'purF', 'purK', 'purL', 'purM~~~purM_2~~~purM_1', 'purQ', 'purR', 'putA~~~putA_2~~~putA_1', 'putP_2~~~putP~~~putP_1', 'pycA', 'pykA', 'pyrAA', 'pyrB_2~~~pyrB', 'pyrC', 'pyrD', 'pyrE', 'pyrF', 'pyrR', 'qoxA', 'qoxB', 'qoxC', 'queC~~~queC_2~~~queC_1', 'queE', 'queH', 'radA', 'rarA', 'rarD~~~rarD_2~~~rarD_1', 'rbgA', 'rbsD', 'rbsK', 'rbsU_1~~~rbsU_2~~~rbsU', 'recD', 'recF', 'recG', 'recN_1~~~recN_2~~~recN', 'recR', 'recU_1~~~recU', 'recX', 'relA', 'rex', 'rhaR~~~araC_2~~~araC_1', 'rhaS~~~mtrA', 'rho', 'ribA', 'ribB', 'ribC_1~~~ribC_2~~~ribC', 'ribH', 'rihA~~~rihA_1~~~rihA_2', 'rimI', 'rimL', 'rimM', 'rimP', 'rlmN_1~~~rlmN_2~~~rlmN', 'rluA1', 'rluB', 'rluD~~~rluD_1~~~rluD_2', 'rnc', 'rnhB', 'rnhC_2~~~rnhC_1~~~rnhC', 'rnr', 'rocA_1~~~rocA_2~~~rocA', 'rpe', 'rpiA_1~~~rpiA_2~~~rpiA', 'rplA', 'rplB', 'rplC_3~~~rplC_2~~~rplC', 'rplD', 'rplF', 'rplJ', 'rplK', 'rplO', 'rplY', 'rpoB', 'rpoC', 'rpoE', 'rpsB', 'rpsC', 'rpsD_1~~~rpsD_2~~~rpsD', 'rpsE', 'rpsG', 'rpsL', 'rsbU_1~~~rsbU_2~~~rsbU', 'rsgA', 'ruvA~~~ruvA_2~~~ruvA_1', 'saeR~~~saeR_1~~~saeR_2', 'saeS~~~sasA~~~saeS_2~~~saeS_1', 'sasC_1~~~sasC_3~~~sasC_2~~~sasC_4~~~sasC', 'sasF', 'sasH~~~sasH_2~~~sasH_1', 'sat_2~~~sat_1~~~sat', 'sceD', 'scoA~~~scoA_2~~~scoA_1', 'scoB~~~scoB_1~~~scoB_2', 'scpA_2~~~scpA_1~~~scpA', 'scpB', 'scrB', 'scrR', 'sdaAA', 'sdaAB', 'sdhA', 'sdhB', 'sdhC_1~~~sdhC_2~~~sdhC', 'sdrH', 'secA', 'secA2_3~~~secA2_2~~~secA2_1~~~secA2', 'secF', 'secY_2', 'secY_3~~~secY_2~~~secY_1', 'sepF', 'serA~~~serA_2~~~serA_1', 'serS', 'sgtB', 'slyA', 'smbA', 'smc_2~~~smc_1~~~smc~~~smc_3', 'smpB', 'speG_1~~~speG_2~~~speG', 'spoUtrmH~~~spoU', 'srrA', 'srtA', 'srtB', 'sua5', 'sucC', 'sucD', 'sufB', 'sufC', 'sufD', 'sufS_1', 'sufS_2~~~sufS_1~~~csd', 'suhB', 'sun', 'tagB', 'tagF', 'tagG', 'tagH', 'tagX', 'tarF~~~tarL', 'tatC_1~~~tatC~~~tatC_2', 'tatP', 'tcaA', 'tcaB_1~~~bcr_2~~~tcaB_2~~~tcaB', 'tcaB_2~~~sotB', 'tdh~~~adh', 'tdk', 'telA_1~~~telA', 'tenA', 'terC', 'tetR_2~~~tetR_3~~~tetR_1', 'tetX_3~~~tetX_2~~~tetX_1', 'tetX~~~tetX_1~~~tetX_2', 'tgt', 'thiD1_1~~~thiD1_2~~~thiD1', 'thiE_2', 'thiG', 'thiI', 'thiM', 'thiN', 'thiO', 'thrA_3~~~thrA_2~~~thrA_1~~~thrA', 'thrB~~~thrB_2~~~thrB_1', 'thrC~~~thrC_2~~~thrC_1', 'thrS', 'thyA_1', 'tig', 'tktA_2~~~tktA~~~tktA_1', 'tmcAL', 'tmk', 'topA~~~topA_1', 'topB_1~~~topB_2~~~topB', 'tpi', 'tpx', 'trkA', 'trkH_2~~~trkH_1~~~trkH', 'trmB_2~~~trmB_1~~~trmB', 'trmD', 'trmE', 'trmU', 'trpA', 'trpB~~~trpB_2~~~trpB_1', 'trpD~~~trpD_1~~~trpC_1~~~trpC_2~~~trpC', 'trpE_1~~~trpE~~~trpE_2', 'trpF', 'trpG_2~~~trpG_1~~~trpG', 'truA', 'truB', 'trxB', 'tsf', 'tuf~~~tuf_2', 'typA', 'tyrA_2~~~tyrA~~~tyrA_1', 'udk', 'ung', 'uppS~~~uppS_1', 'upp_1~~~upp~~~upp_2', 'uraA_1~~~uraA_2~~~uraA', 'ureB', 'ureC_1~~~ureC_2~~~ureC', 'ureE', 'ureF_2~~~ureF_1~~~ureF', 'ureG_1~~~ureG_2~~~ureG', 'uspA_2', 'utp_2~~~utp_1~~~utp', 'uvrA_1', 'uvrA_1~~~uvrA_2', 'uvrB', 'uvrC', 'valS_1~~~valS~~~valS_2', 'vraA_1~~~vraA_2~~~vraA', 'vraB~~~vraB_2~~~vraB_1', 'vraD_2~~~vraD_1~~~vraD', 'vraR', 'vraS_1~~~vraS_2~~~vraS', 'walK', 'walR', 'wecG', 'whiA', 'xerC~~~xerC_2~~~xerC_1', 'xerD', 'xprT', 'yacO', 'yacP', 'yagU', 'yceI_1~~~yceI_2~~~yceI', 'ychF', 'ycnJ', 'ydaD', 'ydjF~~~deoR_1', 'yeeO', 'yehR', 'yfhO_1~~~yfhO_3~~~yfhO_2~~~yfhO_4', 'yfhO_2~~~yfhO_3~~~yfhO_1~~~yfhO', 'yfhP', 'yfiA', 'yflT', 'yicL', 'yidK', 'yjeE', 'ykoD', 'yokF', 'yqfL_1', 'yrrB', 'yteP_2~~~yteP_1~~~yteP', 'yycB_1~~~yycB_2~~~yycB', 'zupT~~~zupT_1~~~zupT_2']

df_group = pd.read_excel(r"D:\python\3.MLST_SC\5.MLST\2.groups\gruops_tree_input.xlsx", index_col=0)
ds = df_group.group

for gene in genes_list:
    file_path = os.path.join(path_input, f"{gene}.aln.fas")
    df = align_df(file_path)
    # df的总长度要大于length
    if len(df.columns) >= 400:
        df_result = best_slide(df, ds, 400, gene)
        csv_ouptut = os.path.join(path_output, f"{gene}.csv")
        df_result.to_csv(csv_ouptut)