import os

ctrl=open('ice-self-similar-aggregates_template.cfg','r').readlines()

freqs=[]

freqs=[13.8,35.5]

for freq in freqs:
    fout=open('control.cfg','w')
    for l in ctrl:
        if 'wavelength' in l:
            l='wavelength %f'%(300/freq/1000.)
        if 'output' in l and '.nc' in l:
            l='output ice-self-similar-aggregates_%06.2f-GHz_scat.nc'%freq
        fout.write(l)
    fout.close()
    os.system('../../src/scatter control.cfg')


ctrl=open('liquid-water_template.cfg','r').readlines()

for freq in freqs:
    fout=open('control%f.cfg'%freq,'w')
    for l in ctrl:
        if 'wavelength' in l:
            l='wavelength %f'%(300/freq/1000.)
        if 'output' in l and '.nc' in l:
            l='output liquid-water_%06.2f-GHz_scat.nc'%freq
        fout.write(l)
    fout.close()
    os.system('../../src/scatter control%f.cfg'%freq)
