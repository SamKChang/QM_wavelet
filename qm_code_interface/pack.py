import qctoolkit as qtk
import numpy as np
import glob

setting = {'program': 'gaussian', 'unit': 'kcal'}
outs_free = qtk.pload('freeAtoms.pkl')
#outs_free = {}
#outs_free['H'] = qtk.QMOut('freeAtoms/H1/H1.out', **setting)
#outs_free['F'] = qtk.QMOut('freeAtoms/F1/F1.out', **setting)

def get_E_data(outs, index=None):
    global outs_free
    E_out = {
        'Ea': [],
        'Ekin': [],
        'Eext': [],
        'Eee': [],
        'Ex': [],
        'Enn': [],
    }
    if index is not None:
        E_out['bond_lengths'] = []
    for o in outs:
        
        # calculate atomic contribution
        stoich = o.molecule.stoichiometry(output='dictionary')
        atoms = {
            'E': 0,
            'Eext': 0,
            'Eee': 0,
            'Ekin': 0,
            'Ex': 0,
        }
        for atom, count in stoich.iteritems():
            atoms['E'] += outs_free[atom].Et * count
            if atom != 'H':
                atoms['Eee'] += outs_free[atom].energies['Eee'] * count
                atoms['Ex'] += outs_free[atom].energies['Ex'] * count
            atoms['Eext'] += outs_free[atom].energies['Eext'] * count
            atoms['Ekin'] += outs_free[atom].energies['Ekin'] * count
            
        # get each energy component    
        E_out['Ea'].append(o.Et - atoms['E'])
        for key in o.energies.keys():
            if key in atoms:
                E_out[key].append(o.energies[key] - atoms[key])
        if 'Enn' in o.energies.keys():
            E_out['Enn'].append(o.energies['Enn'])
            
        if index is not None:
            d = []
            for ind in index:
                i, j = ind
                d.append(o.molecule.distance(i, j))
            E_out['bond_lengths'].append(d)
            
    if index is not None:
        E_out['bond_lengths'] = np.array(E_out['bond_lengths'])

    return E_out

def pack(data_in, E_in):
    data = qtk.ML.pack(data_in)
    for key in E_in.keys():
        if key == 'Ea':
            data['E'] = np.array(E_in[key])
        elif '_inter' not in key:
            data[key] = np.array(E_in[key])
    return data


outs_list = []
for f in sorted(glob.glob('wavelet_HX*_relaxed_udft')):
  outs = []
  for o in sorted(glob.glob('%s' % f + '/*/*.out')):
    out = qtk.QMOut(o, **setting)
    if not np.isnan(out.Et):
      outs.append(out)
  outs_list.append(outs)

names = ['data_HX%d_relaxed.pkl' % m for m in (2, 3,4,5,6)]
for i in range(len(outs_list)):
  outs = outs_list[i]
  print outs[0].name,
  print names[i]
  E = get_E_data(outs)
  data = pack(outs, E)
  qtk.pdump(data, names[i])
