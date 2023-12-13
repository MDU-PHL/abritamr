import pandas, pathlib, numpy, subprocess,logging,json,datetime
from abritamr.CustomLog import CustomFormatter

logger =logging.getLogger(__name__) 
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(CustomFormatter())
fh = logging.FileHandler('update_abritamr_db.log')
fh.setLevel(logging.INFO)
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
fh.setFormatter(formatter)
logger.addHandler(ch) 
logger.addHandler(fh)

def _get_date():

    return datetime.datetime.today().strftime('%Y-%m-%d')

def _get_keys(cfg):

    rename_key = cfg['rename_key']
    other_amr = cfg['other_amr']
    other_non_amr = cfg['other_non_amr']
    oxa_phen_list = cfg['oxa_phen_list']
    adrs = cfg['email_address']
    return rename_key,other_amr,other_non_amr,oxa_phen_list,adrs

def _get_vars():

    cfg = pathlib.Path(__file__).parent / "db" / "update_vars.json"
    if cfg.exists():
        with open(f"{cfg}", "r") as j:
            return _get_keys(json.load(j))
    else:
        logger.critical(f"It seems that {cfg} does not exist. Please check your installation and try again.")
        raise SystemExit
    
def _archive_old_ref_catalog():
    tdy = _get_date()
    og_ref = pathlib.Path(__file__).parent / "db" / "refgenes_latest.csv"
    archive = pathlib.Path(__file__).parent / "db" / f"refgenes_latest.csv.{tdy}"
    logger.info(f"Archiving reference catalog to : {archive}")
    subprocess.run(["cp", f"{og_ref}", f"{archive}"])

def _get_new_catalog():
    logger.info(f"Getting reference catalog from ncbi.")
    updated_html = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"
    p = subprocess.run(f"wget {updated_html}", shell = True, capture_output = True, encoding = "utf-8")
    if p.returncode == 0:
        logger.info(f"Reference catalog downloaded.")
        return True
    else:
        logger.critical(f"Something went wrong with the download of reference catalog. {p.stderr}.")
        raise SystemExit

def _make_key(df):
    if 'key' not in list(df.columns):
        df['mut_acc'] = df[['allele','whitelisted_taxa','refseq_nucleotide_accession']].apply( lambda x: '_'.join(x), axis = 1)
        df['key'] = numpy.where(df['refseq_protein_accession'] == '',df['genbank_protein_accession'],df['refseq_protein_accession'])
        df['key'] = numpy.where(df['subtype'] == 'POINT', df['mut_acc'], df['key'])

    return df

def _open_catalog():
    _get_new_catalog()
    logger.info(f"Generating a dataframe for update.")
    _new = pandas.read_csv("ReferenceGeneCatalog.txt", sep = "\t")
    _new = _new.fillna('')
    _new = _make_key(df = _new)
    return _new

def _get_previous_refgenes():

    existing = _check_existing()

    if existing != False:
        logger.info(f"An existing refgenes file has been found.")
        df = pandas.read_csv(existing)
        df = df.fillna('')
        return pandas.read_csv(existing)
    else:
        logger.warning(f"There is not an existing refgenes file. Refgenes file will be created from new.")
    return existing


def _check_existing():

    existing = pathlib.Path(__file__).parent /"db" / "refgenes_latest.csv"
    
    if existing.exists():
        return existing
    return False

def _capitalise(x):

    a = x.split('/')
    a = [i.capitalize() for i in a]
    _list = list(map(lambda x: x.replace('Carbapenem','Carbapenemase'), a))
 
    return '/'.join(_list)

def cfr(row):

    if row['gene_family'] == 'cfr' and row['class'] == '':
        # print(row)
        return 'Multidrug', _capitalise(row['subclass'])
    else:
        return _capitalise(row['class']),_capitalise(row['subclass'])


def _aminoglycosides(row):

    if 'rRNA' in row['product_name'] and 'methyltransferase' in row['product_name'] and 'AMINOGLYCOSIDE' in row['class']:

        return _capitalise(row['class']),'Aminoglycosides (Ribosomal methyltransferase)'

    else:
        return _capitalise(row['class']),_capitalise(row['subclass'])


def _oxa_phen(row):

    return 'Oxazolidinone/Phenicol',_capitalise(row['subclass'])

def _beta_lactams(row):
    
    if row['subtype'] == 'AMR-SUSCEPTIBLE':
        return 'Beta-lactam','Beta-lactam (Susceptible)'
    elif 'carbapenem-hydrolyzing' in row['product_name'] and 'OXA-51' not in row['product_name'] and row['subtype'] != 'POINT':
        # carbapenem-hydrolyzing
        return 'Beta-lactam','Carbapenemase'
    elif 'carbapenem-hydrolyzing' in row['product_name'] and 'OXA-51' not in row['product_name'] and row['subtype'] == 'POINT':
        # carbapenem-hydrolyzing
        return 'Beta-lactam','Carbapenemase'
    elif 'metallo-beta-lactamase' in row['product_name'] and row['subclass'] == 'CARBAPENEM':
        return 'Beta-lactam','Carbapenemase (MBL)'
    elif 'OXA-51 family' in row['product_name'] or row['allele'] == 'OXA-51':
        return 'Beta-lactam','Carbapenemase (OXA-51 family)'
    elif row['gene_family'] == 'blaZ':
        return 'Beta-lactam','Penicillin resistance (Staphylococcus aureus)'
    elif 'class C' in row['product_name']:
        return 'Beta-lactam','AmpC'
    elif row['gene_family'] == 'blaKPC' and row['subclass'] != 'CARBAPENEM':
        return 'Beta-lactam', 'ESBL (KPC variant)'
    # elif 'OXA-48 family' in row['product_name']:
    #     return 'Beta-lactam',_capitalise(row['subclass'])
    elif not 'class C' in row['product_name'] and ('extended-spectrum' in row['product_name'] or 'extended spectrum' in row['product_name']):
            return 'Beta-lactam','ESBL'
    else:
        return _capitalise(row['class']),_capitalise(row['subclass'])

def _other_antimicrobials(row):

    return 'Other antimicrobial', _capitalise(row['subclass'])

def _other_non_antimicrobials(row):

    return 'Other non-antimicrobial', _capitalise(row['subclass'])

def virulence(row):

    if row['class'] == 'INTIMIN':
        return 'Virulence', f"{_capitalise(row['class'])}_{row['subclass'].lower()}"
    elif 'STX' in row['class']:
        return 'Virulence', _capitalise(row['subclass'])
    elif row['class'] == '':
        return 'Virulence', 'Other'
    else:
        return 'Virulence', _capitalise(row['subclass'])

def _rename(rename_key, row):

    cls = _capitalise(row['class']) if row['class'] not in ['FLUOROQUINOLONE',"NITROFURAN"] else rename_key[row['class']]
    sbcls = rename_key[row['subclass']] if row['subclass'] in rename_key else _capitalise(row['subclass'])

    return cls,sbcls

def _logic(_dict,other_amr,other_non_amr,rename_key):
    logger.info(f"Applying classification logic.")
    for row in _dict:
        if row['class'] == 'AMINOGLYCOSIDE':
            row['enhanced_class'],row['enhanced_subclass'] = _aminoglycosides(row)
        # elif row['class'] == 'PHENICOL/OXAZOLIDINONE' or row['gene_family'] in oxa_phen_list:
        #     row['enhanced_class'],row['enhanced_subclass'] = _oxa_phen(row)
        elif row['gene_family'].startswith('cfr'):
            row['enhanced_class'],row['enhanced_subclass'] = cfr(row)
        elif row['class'] == 'BETA-LACTAM':
            row['enhanced_class'],row['enhanced_subclass'] = _beta_lactams(row)
        elif row['class'] in other_amr:
            row['enhanced_class'],row['enhanced_subclass'] = _other_antimicrobials(row)
        elif row['type'] == 'VIRULENCE':
            row['enhanced_class'],row['enhanced_subclass'] = virulence(row)
        elif row['class'] in other_non_amr:
            row['enhanced_class'],row['enhanced_subclass'] = _other_non_antimicrobials(row)
        elif row['class'] in rename_key:
            row['enhanced_class'],row['enhanced_subclass'] = _rename(rename_key=rename_key, row=row)
        elif row['class'] == 'MULTIDRUG':
            row['enhanced_class'],row['enhanced_subclass'] = 'Multidrug',_capitalise(row['subclass'])
        elif 'multidrug' in row['product_name'] and  (row['class'] == '' and row['subclass'] == ''):
            row['enhanced_class'],row['enhanced_subclass'] = 'Multidrug','Other'
        elif row['type'] != 'VIRULENCE' and (row['class'] == '' and row['subclass'] == ''):
            row['enhanced_class'],row['enhanced_subclass'] = 'Other','Other'
        else:
            row['enhanced_class'],row['enhanced_subclass'] = _capitalise(row['class']),_capitalise(row['subclass'])
    return _dict

def _make_dict(df,other_amr,other_non_amr,rename_key):
    _new_dict = df.to_dict(orient='records')
    _new_dict = _logic(_dict = _new_dict,other_amr =other_amr,other_non_amr=other_non_amr,rename_key=rename_key)
    
    return _new_dict

def _email(adrs,pth):
    # f"cat {_email} | mailx -s '{runid} QC' {a} {user_email}"

    cmd = f"echo 'Please confirm curation of the abritAMR refgenes database attached here.' | mailx -s 'abritamr update' -a {pth} {adrs}"
    
    p = subprocess.run(cmd, shell = True, capture_output=True, encoding='utf-8')

    if p.returncode == 0:
        logger.info("abritAMR update email sent for approval.")
    else:
        logger.critical(f"Something has gone wrong. {p.stderr}.")

def _save_df(df):
    
    sfx = _get_date()
    df.to_csv(f"refgenes_{sfx}.csv", index = False)
    return f"refgenes_{sfx}.csv"

def _updated_entries(new_catalog, previous_catalog):
    # print(previous_catalog)
    previous_catalog = previous_catalog.rename(columns = {"class_new": "class_prev", "subclass_new":"subclass_prev"})
    new_catalog = new_catalog.rename(columns = {"class":"class_new", "subclass": "subclass_new"})
    _tmp= new_catalog.merge(previous_catalog, on = ['key'])
    print(_tmp.columns)
    _tmp['changed'] = numpy.where((_tmp['class_prev'] != _tmp['class_new']) , 'updated', '')
    _tmp['changed'] = numpy.where((_tmp['subclass_prev'] != _tmp['subclass_new']) , 'updated', _tmp['changed'])
    
    _tmp['Previous_class'] = numpy.where(_tmp['changed'] == 'updated', _tmp['class_prev'], '')
    _tmp['Previous_subclass'] = numpy.where(_tmp['changed'] == 'updated', _tmp['subclass_prev'],'')
    changed = list(_tmp[_tmp['changed']!='']['key'])
    
    new_catalog['Status'] = numpy.where(new_catalog['key'].isin(changed), 'updated',new_catalog['Status'])
    new_catalog = new_catalog.merge(_tmp[['key','Previous_class','Previous_subclass']], on = ['key'], how = 'left')
    return new_catalog

def _new_entries(new_catalog, previous_catalog):
    
    logger.info(f"Checking for new entries.")
    new_catalog['Status'] = numpy.where(new_catalog['key'].isin(list(previous_catalog['key'])), 'existing','new')
    logger.info(f"Checking for updated entries.")
    new_catalog = _updated_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)
    return new_catalog

def _update_status(new_catalog, previous_catalog):

    new_catalog = _new_entries(new_catalog=new_catalog,previous_catalog=previous_catalog)

    return new_catalog

def _compare_to_existing(new_catalog, previous_catalog):

    tab = pandas.DataFrame(new_catalog)
    
    if isinstance(previous_catalog,pandas.DataFrame):
        new_catalog = _update_status(new_catalog=tab,previous_catalog=previous_catalog)
        return new_catalog
    else:
        logger.info(f"There is in no previous file to compare with.")
        return tab
    

def create_refgenes():

    rename_key,other_amr,other_non_amr,oxa_phen_list,adrs = _get_vars()
    new_catalog_df = _open_catalog()
    new_catalog = _make_dict(df = new_catalog_df ,other_amr =other_amr,other_non_amr=other_non_amr,rename_key=rename_key)
    existing = _get_previous_refgenes()
    new_catalog = _compare_to_existing(new_catalog=new_catalog,previous_catalog=existing)
    pth = _save_df(df=new_catalog)
    logger.info(f"Saved new refgenes catalog as {pth} and emailing for confirmation.")
    _email(adrs = adrs, pth = pth)
    
