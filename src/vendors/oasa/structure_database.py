#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------

import os, sys, re
from . import oasa_exceptions
try:
    from pysqlite2 import dbapi2 as sqlite
except ImportError as e:
    raise Exception( "The required pysqlite module could not be loaded. More info here: '%s'" % e)

from . import inchi

class Config:
    database_file = os.path.abspath( os.path.join( os.path.dirname( __file__), "structures.db"))


def create_database():
  connection = sqlite.connect( Config.database_file)
  c = connection.cursor()
  c.execute( "DROP TABLE IF EXISTS structures;")
  c.execute( """CREATE TABLE structures (
  id INTEGER PRIMARY KEY,
  name TEXT,
  inchikey TEXT,
  smiles TEXT);""")
  connection.commit()
  c.execute( """CREATE TABLE synonyms (
  id INTEGER,
  synonym TEXT);""")
  connection.commit()
  c.close()
  connection.close()


def normalize_inchi( inchi):
    if inchi.startswith("InChI="):
        return inchi[6:]
    else:
        return inchi

def _open_infile( infilename):
    if os.path.isfile( infilename):
        if infilename.endswith(".gz"):
            import gzip
            f = gzip.open( infilename, "r")
        else:
            f = file( infilename, "r")
        return f
    else:
        return None

def fill_database( infilename, name_cutoff=26, atom_count_cutoff=100):
    """compound is added if either name length of atom_count is below the
    corresponding cutoff value"""
    added = 0
    ignored = 0
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    f = _open_infile( infilename)
    if not f:
        raise ValueError( "File does not exist:", infilename)
    f.readline() # skip the first line
    i = 0
    for line in f:
        if i % 10000 == 0:
            print("done %8d, added %8d, ignored %8d" % (i,added,ignored))
        values = [x.strip() for x in line.strip().split("\t")]
        if len( values) != 4:
            print("Ignoring line:", line, file=sys.stderr)
            continue
        cid, inchikey, smiles, name = [x.strip() for x in line.strip().split("\t")]
        c.execute( "DELETE FROM structures WHERE id=?", (cid,))
        if len( name) <= name_cutoff and len( [x for x in smiles if x.isupper()]) <= atom_count_cutoff and _allow_molecule( name, smiles):
            c.execute( "INSERT INTO structures (id,name,inchikey,smiles) VALUES (?,?,?,?);", (cid, name, inchikey, smiles))
            added += 1
        else:
            ignored += 1
        if i % 100 == 0:
            connection.commit()
        i += 1
    f.close()
    connection.commit()
    connection.close()
    return added, ignored


def get_compounds_from_database( database_file=None, **kw):
    """easy to use interface to the SQL - keyword arguments are converted to
    corresponding SQL commands.
    Examples:
    get_compounds_from_database( smiles='C1CCCCC1')
    get_compounds_from_database( inchi='1/C4H10/c1-3-4-2/h3-4H2,1-2H3')    
    """
    for fname in (database_file, Config.database_file):
        if fname and os.path.exists(fname):
            break
    else:
        raise oasa_exceptions.oasa_error("Structure database not found. Try running 'python structure_database.py structures.txt.gz' in oasa directory to create the database file from default source.")

    if 'inchi' in kw:
        if not 'inchikey' in kw:
            from . import inchi_key
            kw['inchikey'] = inchi_key.key_from_inchi( kw['inchi'])
        del kw['inchi']
    search = ["%s=?" % k for k in kw.keys()]
    values = kw.values()
    tables = "structures"
    if "synonym" in kw:
        tables = "structures, synonyms"
        search.append( "structures.id=synonyms.id")
    if search:
        sql = "SELECT structures.* FROM %s WHERE %s" % (tables, " AND ".join( search))
    else:
        sql = "SELECT structures.* FROM %s" % (tables)
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    try:
        c.execute( sql, values)
    except sqlite.OperationalError as e:
        raise oasa_exceptions.oasa_error( "Error reading from structure database: '%s'" % e)
    ret = []
    c2 = connection.cursor()
    for row in c:
        c2.execute( "SELECT synonym FROM synonyms WHERE id=%d" % row[0])
        synonyms = list( c2)
        ret.append( row+([x[0] for x in synonyms],))
    c.close()
    c2.close()
    connection.close()
    return ret

def find_molecule_in_database( mol, database_file=None):
    """tries to find oasa.molecule mol in the database by using its InChiKey"""
    inchikey = inchi.generate_inchi_key( mol)[0]
    res = get_compounds_from_database( inchikey=inchikey)
    return res

def _allow_molecule( name, smile):
    if smile.count("-]") > 2:
        # more than 2 negative charges
        return False
    if smile.count("+]") > 2:
        # more than 2 positive charges
        return False
    if re.search( r"\[\d", smile):
        return False
    if smile.count(".") > 2:
        return False
    disallowed = ["Tl","U","Ra","Rn","At","Po","Bi","Ta","W","Hf","La","Y","Sr","Nb","Zr","Eu","Lu","Ho","Ce","Nd","Dy","Th","Gd","Yb","Sm",'Es','Np','Pu','Pr','Md','Pm','Pa','Tb','Bk','Cm','Er','Cf','Fm','No',"Am","Ac","Tm","Lr","Fr"]
    for ch in disallowed:
        if ch in smile:
            return False
    # try smiles
##     try:
##         smiles.text_to_mol( smile, calc_coords=False)
##     except oasa_exceptions.oasa_invalid_atom_symbol as e:
##         return False
    return True
    
def filter_src_file( infilename, name_cutoff=26, atom_count_cutoff=1000):
    """print only those lines in infilename that would be allowed
    in structure_database using current settings"""
    added = 0
    ignored = 0
    f = _open_infile( infilename)
    if not f:
        raise ValueError
    f.readline() # skip the first line
    for line in f:
        values = [x.strip() for x in line.strip().split("\t")]
        if len( values) != 4:
            print("Ignoring line:", line, file=sys.stderr)
            continue
        cid, inchikey, smiles, name = [x.strip() for x in line.strip().split("\t")]
        if len( name) <= name_cutoff and len( [x for x in smiles if x.isupper()]) <= atom_count_cutoff and _allow_molecule( name, smiles):
            print(line[:-1])
            added += 1
        else:
            ignored += 1
    f.close()
    return added, ignored

def _is_cid_in_db( cid, cursor):
    cursor.execute( "SELECT id FROM structures WHERE id=%s" % cid)
    if list( cursor):
        return True
    else:
        return False

def add_synonyms_old( fname, only_first=3):
    """this version checks each cid against databse, one by one,
    this is slow, but requires less memory"""
    f = _open_infile( fname)
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    c2 = connection.cursor()
    i = 0
    last_cid = None
    line_count = 0
    for line in f:
        line_count += 1
        parts = line.split()
        cid = parts[0]
        if last_cid != cid:
            last_cid == cid
            count = 1
            exists = _is_cid_in_db( cid, c2)
        else:
            count += 1
        if exists and count <= only_first:
            name = " ".join( parts[1:])
            c.execute( "INSERT INTO synonyms (id,name) VALUES (?,?);", (cid, name))
            i += 1
            if i % 500 == 0:
                connection.commit()
            if i % 10000 == 0:
                print("Added %d" % i)
        if line_count % 100000 == 0:
            print("-- %dk lines --" % (line_count // 1000))
        if line_count % 1000000 == 0:
            break
    connection.commit()
    return i    


def add_synonyms( fname, only_first=3):
    """this version takes all cids into the memory for later comparison and is therefore
    much faster, it can however consume much memory for big databases"""
    f = _open_infile( fname)
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    c2 = connection.cursor()
    i = 0
    last_cid = None
    line_count = 0
    c.execute( "SELECT id FROM structures;")
    cids = set( [row[0] for row in c])
    for line in f:
        line_count += 1
        parts = line.split()
        cid = int( parts[0])
        if last_cid != cid:
            last_cid = cid
            count = 1
            exists = cid in cids
        else:
            count += 1
        if exists and count <= only_first:
            name = " ".join( parts[1:])
            c.execute( "INSERT INTO synonyms (id,synonym) VALUES (?,?);", (cid, name))
            i += 1
            if i % 500 == 0:
                connection.commit()
            if i % 10000 == 0:
                print("Added %d" % i)
        if line_count % 1000000 == 0:
            print("-- %dk lines --" % (line_count // 1000))
        #if line_count % 1000000 == 0:
        #    break
    connection.commit()
    return i    



if __name__ == "__main__":
    import sys
##     a,i = filter_src_file( 'selected_names_name50_inchi100.clean.txt')
##     print >> sys.stderr, a, i
##     sys.exit()

    from optparse import OptionParser
    op = OptionParser(usage="python %prog -c command [options] [structure file]")

    op.add_option( "-c", "--command", action="store",
                   dest="command", default="test",
                   help="command - what to do; one of 'test','update','rebuild','synonyms'")
    op.add_option( "-s", "--synonyms", action="store",
                   dest="synonyms", default="",
                   help="file containing synonyms")
    (options, args) = op.parse_args()

    if options.command == "test":
        print(get_compounds_from_database( inchi="1/C4H10/c1-3-4-2/h3-4H2,1-2H3"))
        print(get_compounds_from_database( smiles="O"))
        from . import smiles
        print(find_molecule_in_database( smiles.text_to_mol("C1CCC=CC1")))
        print(find_molecule_in_database( smiles.text_to_mol("c1ccccc1C")))
        print(get_compounds_from_database( synonym="toluene"))
    elif options.command in ('update','rebuild','synonyms'):
        if not options.command == 'synonyms':
            if len( args) >= 1:
                fname = args[0]
            else:
                print("You must supply a valid filename of a file containing structures to be read.", file=sys.stderr)
                sys.exit()
            if options.command == 'rebuild':
                if os.path.exists( Config.database_file):
                    try:
                        os.remove( Config.database_file)
                    except Exception as e:
                        print("Old database file could not be removed. Reason: %s" % e, file=sys.stderr)
                        print("Quitting..", file=sys.stderr)
                        sys.exit()
                create_database()
            added, ignored = fill_database( fname)
            print("Added: %d, Ignored: %d" % (added, ignored))
        if options.synonyms:
            print("Adding synonyms")
            #import profile
            #profile.run( 'add_synonyms( options.synonyms)')
            i = add_synonyms( options.synonyms)
            #i = 0
            print("Added %d synonyms" % i)
    else:
        if options.command not in ("test","rebuild","update"):
            print("Invalid action '%s'" % options.command, file=sys.argv)
