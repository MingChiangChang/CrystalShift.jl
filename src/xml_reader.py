from pathlib import Path

import xml.etree.ElementTree as ET
import numpy as np
import re
home = Path.home()

path = home / 'Desktop' / 'AlLiFe_data' / 'sticks' / 'AlLiFe_oxides'

xmls = list(path.glob("*.xml"))

with open('sticks.csv', 'w') as f:

    for idx, xml in enumerate(xmls):
        tree = ET.parse(xml)
        root = tree.getroot()
        pdf = root.find('pdf_data')
        chem_form = pdf.find('chemical_formula').text
        xstal_sys = pdf.find('xstal_system').text
        sg = pdf.find('xtlsg').text
        print(sg)
        if ' ' in sg:
            sg = sg[:sg.index(' ')]

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
        intensities = []
        for stick in stick_series:
            h = stick.find('h').text
            k = stick.find('k').text
            l = stick.find('l').text
            intensity = stick.find('intensity').text
            intensity = re.sub("[^0-9]", "", intensity)
            hs.append(h)
            ks.append(k)
            ls.append(l)
            intensities.append(intensity)


        intensities = np.array(intensities).astype('float64')
        intensities = intensities/np.max(intensities)
        intensities = intensities * 100

        f.write(f'{idx},{chem_form}_{sg},{xstal_sys},{a},{b},{c},{alpha},{beta},{gamma}')

        for h, k, l, intensity in zip(hs, ks, ls, intensities):
            f.write(f'\n{h},{k},{l},0.0,{intensity}')
        f.write('#\n')