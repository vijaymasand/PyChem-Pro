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



cdml_to_smiles = {'Me': 'C',
                  'Ph': 'c1ccccc1',
                  'Et': 'CC',
                  'OCH3': 'OC',
                  'NO2': 'N(=O)=O',
                  'COOH': 'C(=O)O',
                  'COOCH3': 'C(=O)OC',
                  'CN': 'C#N',
                  'SO3H': 'S(=O)(=O)O',
                  'PPh3': 'P(c1ccccc1)(c2ccccc2)c3ccccc3',
                  'OMe': 'OC',
                  'COCl': 'C(=O)Cl',
                  'CH2OH': 'CO'
                  }


name_to_smiles = {'Me': 'C',
                  'Ph': 'c1ccccc1',
                  'Et': 'CC',
                  'Pr': 'CCC',
                  'Bn': 'Cc1ccccc1',
                  'Bz': 'C(=O)c1ccccc1',
                  'Ac': 'C(=O)C',
                  'Ms': 'S(=O)(=O)C',
                  'Tf': 'S(=O)(=O)C(F)(F)F',
                  'Ts': 'S(=O)(=O)c1ccc(C)cc1',
                  'Bu': 'CCCC',

                  'iPr':   'C(C)C',
                  'tBu':   'C(C)(C)C',
                  'Boc':   'C(=O)OC(C)(C)C',
                  'Alloc': 'C(=O)OCC=C',
                  'Z':     'C(=O)OCc1ccccc1',
                  'Fmoc':  'C(=O)OCC1c2ccccc2-c3c1cccc3',
#                  'DMAM':  '
                  'CE':    'CCC#N',
                  'MOP':   'CC1=CC(OC)=CC=N1=O',
                  'Tr':    'C(c1ccccc1)(c2ccccc2)c3ccccc3',
                  'MMTr':  'C(c1ccc(OC)cc1)(c2ccccc2)c3ccccc3',
                  'DMTr':  'C(c1ccc(OC)cc1)(c2ccc(OC)cc2)c3ccccc3',
                  'TMS':   'Si(C)(C)C',
                  'TBDMS': 'Si(C)(C)C(C)(C)C',
                  'TBDPS': 'Si(C(C)(C)C)(c1ccccc1)c2ccccc2',
                  
                  }
