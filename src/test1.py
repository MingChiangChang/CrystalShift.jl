import numpy
import xrayutilities as xu

if __name__ == '__main__':
    #freeze_support()
    
    path = "/Users/mingchiang/Desktop/github/Crystallography_based_shifting/data/Delta.cif"
    sg = xu.materials.cif.CIFFile(path).SGLattice()
    delta = xu.materials.material.Crystal('Delta', sg)
    powder = xu.simpack.PowderDiffraction(delta)
    # pm = xu.simpack.PowderModel(powder, I0=1)
    print(powder.data.keys())

    #tt = numpy.arange(5, 120, 0.01)
    #Fe_powder = xu.simpack.Powder(xu.materials.Fe, 1,
    #                              crystallite_size_gauss=100e-9)
    #Co_powder = xu.simpack.Powder(xu.materials.Co, 5,  # 5 times more Co
    #                              crystallite_size_gauss=200e-9)
    #pm = xu.simpack.PowderModel(Fe_powder, Co_powder, I0=100)
    #inte = pm.simulate(tt)
    #print(inte)

