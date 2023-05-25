#!/usr/bin/env python
"""
Helper routine to generate a Fortran code snipped that computes the OH
reactivity based on the reactions listed in gckpp_Monitor.F90.
The OH reactivity is defined as the inverse of its lifetime and is
determined using the reaction rates of all reactions that consume OH.
The generated code snipped can then be inserted into one of the Fortran
modules, e.g. gckpp_Util.F90.

Usage:
python OHreact_parser.py

Revision History:
2018-06-18 - christoph.a.keller@nasa.gov - initial version
See git history for subsequent updates
"""

import sys
import datetime


def writeHeader(fo, mechName):
    """
    Writes the header of the Get_OHreactivity subroutine.
    
    Args:
    -----
    fo (_io.TextIOWrapper) : File handle for write command
    """
    fo.write('! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    fo.write('\n')
    fo.write('! Get_OHreactivity - returns the OH reactivity\n')
    fo.write('! The OH reactivity is defined as the inverse of its lifetime.\n')
    fo.write('! This routine was auto-generated using script OHreact_parser.py.\n')
    fo.write('! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    fo.write('\n')
    fo.write('SUBROUTINE Get_OHreactivity ( CC, RR, OHreact )\n')
    fo.write('\n')
    fo.write('! CC - Concentrations of species (local)\n')
    fo.write('  REAL(kind=dp) :: CC(NSPEC)\n')
    fo.write('! RR - reaction rates (local)\n')
    fo.write('  REAL(kind=dp) :: RR(NREACT)\n')
    fo.write('! OHreact - OH reactivity [s-1]\n')
    fo.write('  REAL(kind=dp) :: OHreact\n')
    fo.write('\n')
    if "carbon" in mechName:
        fo.write('  OHreact = 0.0_dp') # carbon mech has no OH reactivity
    else:
        fo.write('  OHreact = ')


def writeFooter(fo):
    """
    Writes the footer of the Get_OHreactivity subroutine.
    
    Args:
    -----
    fo       (_io.TextIOWrapper) : File handle for write command
    """
    fo.write('\n\nEND SUBROUTINE Get_OHreactivity\n')
    fo.write('! End of Get_OHreactivity subroutine\n')
    fo.write('! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')


def writeReactivity(mechName):
    """
    Writes the Get_OHreactivity routine in the gckpp_Util.F90 module.
    
    Args:
    -----
    mechName (str) : Name of the chemical mechanism
    """

    # Filenames
    file_reactions = 'gckpp_Monitor.F90'
    file_species   = 'gckpp_Parameters.F90'
    file_updated   = 'gckpp_Util.F90'

    # Rread file with reactions
    with open(file_reactions, 'r') as f:
        lines = f.readlines()

    # Read files with species
    with open(file_species, 'r') as f:
        params = f.readlines()

    # Read gckpp_Util.F90 and write back out, minus END MODULE line. 
    # Save # that line to write out later.
    EndModuleStr = ''
    with open(file_updated, 'r') as f:
        utillines = f.readlines()
    with open(file_updated, 'w') as f:
        for line in utillines:
            if 'END MODULE' in line:
                EndModuleStr = line
                break
            else:
                f.write(line.rstrip()+'\n')

    # Open gckpp_Util.F90
    fo = open(file_updated,"a")

    # Write header to the Get_OHreactivity routine
    writeHeader(fo, mechName)
    
    nn   = 0
    ntot = 0

    # get all lines with a reaction
    rxt = [i for i in lines if '-->' in i]

    # walk through all reactions
    irct = 0
    for irxt in rxt:
        irct += 1
        irxt = irxt.replace("\n", "")
        irxt = irxt.replace("'", "")
        #---------------------------------------------------
        # Also remove extra text for safety's sake
        #  -- Bob Yantosca (03 May 2022)
        irxt = irxt.replace("/)", "")
        irxt = irxt.replace(",", "")
        irxt = irxt.replace("&", "")
        irxt = irxt.replace("!", "")
        irxt = irxt.replace("up to ", "")
        #---------------------------------------------------
        spl = irxt.split('-->')

        # get reaction number and reactant if OH is on left-hand side of reaction
        if ' OH ' in spl[0]:
            # Go to next line after writing 5 reaction terms
            if nn == 5:
                fo.write(' &\n')
                nn = 0
            # reaction number
            if 'index ' in spl[1]:
                rn = int(spl[1].split('index ')[1])
                if irct != rn:
                    print('Warning: reaction number mismatch')
            else:
                rn = irct
            # species name
            if ' + ' in spl[0]:
                # there is one OH
                nOH = 1
                newspl = spl[0].split('+')
                if ' OH ' in newspl[0]:
                    spc = newspl[1]
                else:
                    spc = newspl[0]
            else:
                if '2 OH' in spl[0]:
                    nOH = 2
                    spc = 'NotAvail'
                else:
                    print('unexpected entry: '+spl)
                    continue
            spc = spc.replace("'","")
            spc = spc.strip()
            # get species index
            if nOH == 1:
                istr = 'ind_'+spc+' ='
                iln = [i for i in params if istr in i]
                if len(iln) != 1:
                    print('cannot match species '+istr)
                    continue
                else:
                    iln = iln[0].replace("\n","")
                    id = int(iln.split('=')[1])

            # now we have the reaction number, species id and number of OH.
            # Can construct string entry
            if nOH == 2:
                istr = '2*RR('+str(rn)+')'
            else:
                istr = 'RR('+str(rn)+')*CC('+str(id)+')'
            if ntot > 0:
                istr = ' + '+istr
                if nn == 0 :
                    istr = '         '+istr

            # write to file
            nn   += 1
            ntot += 1
            fo.write(istr)

    # Close the subroutine and the module
    writeFooter(fo)
    fo.write(EndModuleStr)

    # all done
    print('Reactivity consists of '+str(ntot)+' reactions')
    print('Written to '+file_updated)
    fo.close()


def main():
    """
    Main program.  Gets number of arguments and calls writeReactivity
    """

    # Get the name of the mechansim
    nArgs = len( sys.argv )
    mechName = ""
    if ( nArgs == 2 ): 
        mechName = sys.argv[1]

    # Write the OH reactivity routine in gckpp_Util.F90
    writeReactivity(mechName)
    

if __name__ == '__main__':
    main()
