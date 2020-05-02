""" module to create an ensight compatible file
to visualize your  data"""

import os
import h5py
import numpy as np
from lxml import etree

NSMAP = {"xi": "http://www.w3.org/2001/XInclude"}


# pylint: disable=c-extension-no-member
class NpArray2Xmf():
    """ main class for data output in XDMF format """

    def __init__(self,
                 filename,
                 domain_name=None,
                 mesh_name=None,
                 time=None,
                 xmf_only=False):
        """ class startup"""
        extension = os.path.splitext(filename)
        if extension[-1] == ".xmf":
            self.filename = extension[0] + ".h5"
        elif extension[-1] == ".h5":
            self.filename = filename
        else:
            raise RuntimeError("Only extensions .xmf or .h5 are allowed")

        self.geotype = None
        self.mesh = {}
        self.data = {}
        self.shape = None
        self.mesh["domain"] = domain_name
        self.mesh["mesh"] = mesh_name
        self.mesh["time"] = time
        if time is None:
            self.mesh["time"] = 0.0
        self.mesh["time"] = "%14.8e" % self.mesh["time"]

        if self.mesh["mesh"] is None:
            self.mesh["mesh"] = "Mesh"

        if self.mesh["domain"] is None:
            self.mesh["domain"] = "Domain"

        self.xmf_only = xmf_only

    def create_grid(self, nparray_x, nparray_y, nparray_z):
        """ create the grid according to numpy arrays x, y ,z
        if arrays are 1D, switch to cloud point
        if arrays are 2D, switch to quad connectivity
        if arrays are 3D, switch to hexaedrons connectivity"""
        self.mesh["x"] = np.ravel(nparray_x)
        self.mesh["y"] = np.ravel(nparray_y)
        self.mesh["z"] = np.ravel(nparray_z)
        self.shape = list(nparray_x.shape)
        dim = len(self.shape)
        if dim == 1:
            self.geotype = "cloud"
        if dim == 2:
            self.geotype = "quads"
        if dim == 3:
            self.geotype = "hexas"

        if self.geotype is None:
            raise RuntimeError("Unexpected shape of nparray :"
                               + " ".join(self.shape))

    def add_field(self, nparray_field, variable_name):
        """ add a field, assuming same shape as nparray of coordiantes """
        self.data[variable_name] = nparray_field

    def _type(self, var):
        """ retrun the xmf type according to nparray"""
        var_shape = list(self.data[var].shape)

        dtype = self.data[var].dtype
        numbertype = None
        if dtype == "float64":
            numbertype = "Float"
        if dtype == "S4":
            numbertype = "Char"
        if numbertype is None:
            raise RuntimeError("Array of type " + dtype
                               + "(" + var + ")"
                               + "not recognized")

        attributetype = None
        if var_shape == self.shape:
            attributetype = "Scalar"
        if var_shape[:-1] == self.shape:
            if var_shape[-1] == 3:
                attributetype = "Vector"
        if attributetype is None:
            raise RuntimeError("Var " + var
                               + " of shape " + str(var_shape)
                               + " not consistent  with grid of shape "
                               + str(self.shape) +
                               "\n (neither scalar nor 3D vector...)")

        shape_str = " ".join(str(dim) for dim in self.data[var].shape)
        return (numbertype,
                attributetype,
                shape_str)

    def xmf_dump(self):
        """ create XDMF descriptor file """
        if self.geotype == "cloud":
            topology = "PolyVertex"
        if self.geotype == "quads":
            topology = "2DSMesh"
        if self.geotype == "hexas":
            topology = "3DSMesh"

        dims = " ".join(str(dim) for dim in self.shape)

        xmf_tree = dict()
        xmf_tree['root'] = etree.Element("Xdmf", Version="2.0", nsmap=NSMAP)

        xmf_tree['dom'] = etree.SubElement(xmf_tree['root'], "Domain",
                                           Name=self.mesh['domain'])

        xmf_tree['grd'] = etree.SubElement(xmf_tree['dom'], "Grid",
                                           Name=self.mesh['mesh'],
                                           Type='Uniform')

        etree.SubElement(xmf_tree['grd'], "Time",
                         Type='Single',
                         Value=self.mesh['time'])

        etree.SubElement(xmf_tree['grd'], "Topology",
                         Name="Topo",
                         TopologyType=topology,
                         NumberOfElements=dims)

        xmf_tree['geo'] = etree.SubElement(xmf_tree['grd'], "Geometry",
                                           GeometryType="X_Y_Z")

        for var in ['x', 'y', 'z']:
            field = etree.SubElement(xmf_tree['geo'], "DataItem",
                                     Dimensions=dims,
                                     Format="HDF",
                                     NumberType="Float",
                                     Precision="8")
            text = "%s:/mesh/%s" % (os.path.basename(self.filename), var)
            field.text = "\n%s%s\n%s" % (11 * " ", text, 8 * " ")

        for var in self.data:
            numbertype, attributetype, dims = self._type(var)

            attr = etree.SubElement(xmf_tree['grd'], "Attribute",
                                    Name=var,
                                    Center="Node",
                                    AttributeType=attributetype)

            field = etree.SubElement(attr, "DataItem",
                                     Dimensions=dims,
                                     Format="HDF",
                                     NumberType=numbertype,
                                     Precision="8")
            text = "%s:/variables/%s" % (os.path.basename(self.filename), var)
            field.text = "\n%s%s\n%s" % (11 * " ", text, 8 * " ")

        xmf_file = self.filename.replace(".h5", ".xmf")
        xmf_ct = etree.tostring(xmf_tree['root'],
                                pretty_print=True,
                                xml_declaration=True,
                                doctype='<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
        xmf_ct = xmf_ct.decode().replace("encoding=\'ASCII\'", "")
        with open(xmf_file, "w") as fout:
            fout.write(xmf_ct)

    def dump(self):
        """ dump the final file """
        if not self.xmf_only:
            fout = h5py.File(self.filename, "w")

            mesh_gp = fout.create_group("mesh")
            for coord in ["x", "y", "z"]:
                mesh_gp.create_dataset(coord, data=self.mesh[coord])

            var_gp = fout.create_group("variables")
            for var in self.data:
                var_gp.create_dataset(var, data=self.data[var])
            fout.close()
        self.xmf_dump()


def create_time_collection_xmf(collection_filenames, xmf_filename):
    """ Creates xmf file holding time collection of xmf files

    Parameters :
    ============
    collection_filenames: a list of single time xmf filenames to collect
    xmf_filename : the name of the output file

    Returns:
    ========
    None
    """

    root = etree.Element("Xdmf", Version="2.0", nsmap=NSMAP)
    dom = etree.SubElement(root, "Domain")

    grid = etree.SubElement(dom, "Grid",
                            Name=os.path.split(xmf_filename)[-1],
                            GridType="Collection",
                            CollectionType="Temporal")

    for filename in collection_filenames:
        etree.SubElement(grid, "XI_INCLUDE", href=filename,
                         xpointer='xpointer(//Xdmf/Domain/Grid)')

    xmf_ct = etree.tostring(root,
                            pretty_print=True,
                            xml_declaration=True,
                            doctype='<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    xmf_ct = xmf_ct.decode()
    with open(xmf_filename, "w") as fout:
        xmf_ct = xmf_ct.replace("encoding=\'ASCII\'", "")
        xmf_ct = xmf_ct.replace("XI_INCLUDE", "xi:include")
        fout.write(xmf_ct)


if __name__ == '__main__':
    DIM_X = 41
    DIM_Y = 21
    DIM_Z = 11

    SIZE_X = 4.
    SIZE_Y = 2.
    SIZE_Z = 1.

    # 1D
    TEST_X = np.linspace(0, SIZE_X, DIM_X)
    TEST_Y = np.linspace(0, SIZE_Y, DIM_X)
    TEST_Z = np.linspace(0, SIZE_Z, DIM_X)
    TEST_U = (np.sin(TEST_X / SIZE_X * 1 * np.pi)
              * np.sin(TEST_Y / SIZE_Y * 1 * np.pi)
              * np.sin(TEST_Z / SIZE_Z * 1 * np.pi))

    TEST_F = NpArray2Xmf("./test1D.h5")
    TEST_F.create_grid(TEST_X, TEST_Y, TEST_Z)
    TEST_F.add_field(TEST_U, "foobar")

    TEST_V = np.stack((TEST_U,
                       TEST_U,
                       TEST_U),
                      axis=1)
    TEST_F.add_field(TEST_V, "foobar_vect")

    TEST_F.dump()

    # 2D
    TEST_X = np.tile(np.linspace(0., SIZE_X, DIM_X), (DIM_Y, 1))
    TEST_Y = np.tile(np.linspace(0., SIZE_Y, DIM_Y), (DIM_X, 1)).transpose()
    TEST_Z = np.ones((DIM_Y, DIM_X))
    TEST_U = (np.sin(TEST_X / SIZE_X * 1 * np.pi)
              * np.sin(TEST_Y / SIZE_Y * 1 * np.pi)
              * np.sin(TEST_Z * 0.5 * np.pi))

    TEST_F = NpArray2Xmf("./test2D.h5")
    TEST_F.create_grid(TEST_X, TEST_Y, TEST_Z)
    TEST_F.add_field(TEST_U, "foobar")
    TEST_F.dump()

    TEST_X = TEST_X[:, :, None].repeat(DIM_Z, 2)
    TEST_Y = TEST_Y[:, :, None].repeat(DIM_Z, 2)
    TEST_Z = np.tile(np.linspace(0., SIZE_Z, DIM_Z), (DIM_X, 1)).transpose()
    TEST_Z = TEST_Z[:, :, None].repeat(DIM_Y, 2).transpose((2, 1, 0))
    TEST_U = (np.sin(TEST_X / SIZE_X * 1 * np.pi)
              * np.sin(TEST_Y / SIZE_Y * 1 * np.pi)
              * np.sin(TEST_Z / SIZE_Z * 1 * np.pi))

    TEST_F = NpArray2Xmf("./test3D.h5")
    TEST_F.create_grid(TEST_X, TEST_Y, TEST_Z)
    TEST_F.add_field(TEST_U, "foobar")
    TEST_F.dump()
