import os
import sys
import argparse
import pymzml


protonMass = 1.007276466771

def convert_mz_2_mass( mz, charge ):
    '''
    NOTE: equal to [ charge * mz - ( charge * PROTON) ]
    '''
    return float(charge) * float( float(mz) - float(protonMass) )

def _determine_mzml_name_base( file_name ):
    file_name = os.path.basename( file_name )
    if file_name.upper().endswith('.MZML.GZ'):
        mzml_name_base = file_name[:-8]
    elif file_name.upper().endswith('.MZML'):
        mzml_name_base = file_name[:-5]
    else:
        raise Exception("Can not determine mzml base name from {0}".format(
            file_name
        ))
    return mzml_name_base

def main( mzmlFilename = None, i_decimals = 0, mz_decimals = 4, name_base = '', output_dir = './', mass_decimals = 6):
    if(name_base == ''):
        name_base = _determine_mzml_name_base( mzmlFilename )

    print('Converting file:\n\tmzml : {0}\n\toutput dir : {4}\n\tmass_decimals : {1}\n\tmz_decimals : {2}\n\ti_decimals : {3}\n\tname_base : {5}'.format(
        mzmlFilename,
        mass_decimals,
        mz_decimals,
        i_decimals,
        output_dir,
        name_base
    ))
    
    ## The output directory might not already exist
    if not os.path.exists(output_dir):
        ## Create the output directory
        ## (TODO: should consider using 'mkdirs' instead, if more than 1 level of directory needs to be created)
        os.mkdir(output_dir, 0775)

    ''' Going from mzML to dta: (<ws> represents white space)
           Each "ms level"=2 scan in the mzML file is exported into its own .dta file
           The filename is the original mzML filename without the mzML extension, plus ".<scan start>.<scan end>.<peptide charge state>.dta"
           The first line of the file is: (precursor's 'selected ion m/z')*(charge state) - (charge state - 1)*protonMass <whitespace> <charge state>
           Each subsequent line contains: <m/z value to 4 decimal places> <whitespace> <intensity rounded to nearest integer>
           (Those decimal accuracies are based on what MzXML2Search produces from mzML / mzXML files)
    '''
    '''msrun = pymzml.run.Reader("e007700/e007700_msconvert.mzML")'''
    msrun = pymzml.run.Reader(mzmlFilename)
    numLevel1DTAs = 0
    numLevel2DTAs = 0
    numWrittenDTAs = 0
    numSkippedMissingChargeState = 0
    numSkippedSingleChargeState = 0
    for currSpectrum in msrun:
        'Note: We only create .dta files for ms level 2 scans'
        if currSpectrum['ms level'] != 2:
            numLevel1DTAs += 1
            continue

        scanID = currSpectrum['id']
        extendedScanID = str(scanID).rjust(5, '0')
        numLevel2DTAs += 1;
        if numLevel2DTAs % 500 == 0:
            print(
                'File : {0:^40} : Processing spectrum {1}'.format(
                    os.path.basename( mzmlFilename ),
                    numLevel2DTAs
                )
            )

        m_over_z = float(currSpectrum['precursors'][0]['mz'])
        chargeStateReference = currSpectrum['precursors'][0]['charge']
        if chargeStateReference is None:
            numSkippedMissingChargeState += 1
            continue
        chargeState = int(chargeStateReference)
        if chargeState == 1:
            numSkippedSingleChargeState += 1
            continue

        ''' If we successfully found both chargeState and m/z, lets spit out the peaks to a file '''
        theMass = float(convert_mz_2_mass( m_over_z, chargeState )) + float(protonMass)
        '''
        theMass = (m_over_z)*(chargeState) - (chargeState - 1)*protonMass
        print '%.6f' % theMass
        print 'm_over_z: %f, chargeState: %d' % (m_over_z, chargeState)
        '''
        
        file_name = '{0}.{1}.{1}.{2}.dta'.format(
            name_base,
            extendedScanID,
            chargeState,
            )
        
        filePath = os.path.join(output_dir, file_name)

        '''
        print 'Attempting to write to %s' % filePath
        '''
        
        with open(filePath, 'w') as outfile:
            print >> outfile, '%.*f %d' % ( mass_decimals, theMass, chargeState )

            for mz, intensity in currSpectrum.peaks:
                print >> outfile, '%.*f %.*f' % (
                        mz_decimals,
                        mz,
                        i_decimals,
                        intensity,
                    )
        numWrittenDTAs += 1

    print('Found {0} level 2 scan entries'.format(numLevel2DTAs))
    print('Skipped {0} entries with missing charge state'.format(numSkippedMissingChargeState))
    print('Skipped {0} entries with single charge state'.format(numSkippedSingleChargeState))
    print('Wrote {0} dta entries'.format(numWrittenDTAs))

    
if __name__ == '__main__':
    # parsing command line arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'mzmlFilename',
        help='mzml input file')
    parser.add_argument(
        '-m', '--mass_decimals', default=6,
        help='Number of decimals for MH+ value (default: 6)', type=int)
    parser.add_argument(
        '-mz', '--mz_decimals', default=4,
        help='Number of decimals for m/z values (default: 4)', type=int)
    parser.add_argument(
        '-i', '--i_decimals', default=0,
        help='Number of decimals for intensity values (default: 0)', type=int)
    parser.add_argument(
        '-n', '--name_base', default='',
        help='base name for dta output files (default: original file name)')
    parser.add_argument(
        '-o', '--output_dir', default='./',
        help='output directory for dta files (will be created if doesn\'t exist, default: curr dir.)')

    if len( sys.argv ) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        tmp = main( **args.__dict__ )

