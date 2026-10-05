from pathlib import Path

import xml.etree.ElementTree as ET
import numpy as np
import re
home = Path.home()

path = home / 'Desktop' / 'AlLiFe_data' / 'sticks' / 'AlLiFe_oxides'
path = home / 'Downloads' / 'AlLiFeO_assembled_icdd'
path = home / "Downloads" / "SnO"
#path = home / 'Downloads' / 'AlLiFeO'
#path = home / 'Downloads' / 'CrFeV_toCornell' / 'icdd'
#path = home / 'Downloads' / 'AlLiFeO copy'

xmls = sorted(list(path.glob("*.xml")))

def check_none(*arg):
    for i in arg:
        if i is None:
            return True
    return False

with open(f'{str(path)}/_sticks.csv', 'w') as f:

    for idx, xml in enumerate(xmls):
        print(xml)
        tree = ET.parse(xml)
        root = tree.getroot()
        pdf = root.find('pdf_data')
        chem_form = pdf.find('chemical_formula').text
        xstal_sys = pdf.find('xstal_system').text
        if pdf.find('xtlsg') is not None:
            sg = pdf.find('xtlsg').text
            print(idx, chem_form, sg)
            if ' ' in sg:
                sg = sg[:sg.index(' ')]
        else:
            sg = 'P1'

        a = pdf.find('xtla').text
        b = pdf.find('xtlb').text
        c = pdf.find('xtlc').text
        alpha = pdf.find('xtlal').text
        beta  = pdf.find('xtlbe').text
        gamma = pdf.find('xtlga').text
 
        graph = root.find('graphs')
        stick_series = graph.find('stick_series')

        hs = []
        ks = []
        ls = []
        qs = []
        intensities = []
        for stick in stick_series:
            h = stick.find('h').text
            k = stick.find('k').text
            l = stick.find('l').text
            if check_none(h, k, l):
                continue
            intensity = stick.find('intensity').text
            theta = stick.find('theta').text
            q = 4*np.pi*np.sin(np.radians(float(theta)/2)) / 1.5418
            intensity = re.sub("[^0-9]", "", intensity)
            hs.append(h)
            ks.append(k)
            ls.append(l)
            print(q)
            qs.append(q*10)
            intensities.append(intensity)


        intensities = np.array(intensities).astype(np.float64)
        intensities = intensities/np.max(intensities)
        intensities = intensities * 100

        f.write(f'{idx},{chem_form}_{sg},{xstal_sys},{a},{b},{c},{alpha},{beta},{gamma}')

        for h, k, l, q, intensity in zip(hs, ks, ls, qs, intensities):
            f.write(f'\n{h},{k},{l},{q},{intensity}')
        f.write('#\n')
