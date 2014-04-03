# to sa przyblizone wartosci dla kilku prostych przykladow
# C_x_cor podawana jest dla ix = -1/2, 1/2..., iy = 0, 1..
# C_y_cor podawana jest dla iy = -1/2, 1/2..., ix - 0, 1...

import numpy

wyniki = {}

wyniki["Cx"] = []
wyniki["Cy"] = []
wyniki["psi_0"] = []
wyniki["psi_iter1"] = []
wyniki["psi_iter2"] = []
wyniki["psi_iter3"] = []
wyniki["Ax_iter2"] = []
wyniki["Ay_iter2"] = []
wyniki["Ax_iter3"] = []
wyniki["Ay_iter3"] = []
wyniki["Bx_iter2"] = []
wyniki["By_iter2"] = []
wyniki["Bx_iter3"] = []
wyniki["By_iter3"] = []
wyniki["Cx_cor_iter2"] = []
wyniki["Cy_cor_iter2"] = []
wyniki["Cx_cor_iter3"] = []
wyniki["Cy_cor_iter3"] = []

wyniki["Cx"].append(0.)
wyniki["Cy"].append(1.)
wyniki["psi_0"].append(numpy.array([[0.,0.,0.,],[0.,1.,0.],[0.,0.,0.,]]))
wyniki["psi_iter1"].append(numpy.array([[0.,0.,0.,],[0.,0.,1.],[0.,0.,0.,]]))
wyniki["psi_iter2"].append([wyniki["psi_iter1"][-1]])
wyniki["psi_iter3"].append([wyniki["psi_iter1"][-1]])
wyniki["Ax_iter2"].append(numpy.array([[0.,0.,1.,],[0.,0.,-1.],[0.,0.,0.,]]))
wyniki["Ay_iter2"].append(numpy.array([[0.,0.,0.,],[0.,1.,-1.],[0.,0.,0.,]]))
wyniki["Ax_iter3"].append(wyniki["Ax_iter2"][-1])
wyniki["Ay_iter3"].append(wyniki["Ay_iter2"][-1])
wyniki["Bx_iter2"].append(0.5*numpy.array([[-1.,1.,0.,],[-1.,1.,0.],[0.,0.,0.,]]))
wyniki["By_iter2"].append(0.5*numpy.array([[0.,1.,1.,],[0.,0.,0.],[0.,-1.,-1.,]]))
wyniki["Bx_iter3"].append(wyniki["Bx_iter2"][-1])
wyniki["By_iter3"].append(wyniki["By_iter2"][-1])
wyniki["Cx_cor_iter2"].append(numpy.zeros((3,3)))
wyniki["Cy_cor_iter2"].append(numpy.zeros((3,3)))
wyniki["Cx_cor_iter3"].append(numpy.zeros((3,3)))
wyniki["Cy_cor_iter3"].append(numpy.zeros((3,3)))

wyniki["Cx"].append(1.)
wyniki["Cy"].append(0.)
wyniki["psi_0"].append(numpy.array([[0.,0.,0.,],[0.,1.,0.],[0.,0.,0.,]]))
wyniki["psi_iter1"].append(numpy.array([[0.,0.,0.,],[0.,0.,0.],[0.,1.,0.,]]))
wyniki["psi_iter2"].append([wyniki["psi_iter1"][-1]])
wyniki["psi_iter3"].append([wyniki["psi_iter1"][-1]])
wyniki["Ax_iter2"].append(numpy.array([[0.,0.,0.,],[0.,1.,0.],[0.,-1.,0.,]]))
wyniki["Ay_iter2"].append(numpy.array([[0.,0.,0.,],[0.,0.,0.],[1.,-1.,0.,]]))
wyniki["Ax_iter3"].append(wyniki["Ax_iter2"][-1])
wyniki["Ay_iter3"].append(wyniki["Ay_iter2"][-1])
wyniki["Bx_iter2"].append(0.5*numpy.array([[0.,0.,0.,],[1.,0.,-1.],[1.,0.,-1.,]]))
wyniki["By_iter2"].append(0.5*numpy.array([[-1.,-1.,0.,],[1.,1.,0.],[0.,0.,0.,]]))
wyniki["Bx_iter3"].append(wyniki["Bx_iter2"][-1])
wyniki["By_iter3"].append(wyniki["By_iter2"][-1])
wyniki["Cx_cor_iter2"].append(numpy.zeros((3,3)))
wyniki["Cy_cor_iter2"].append(numpy.zeros((3,3)))
wyniki["Cx_cor_iter3"].append(numpy.zeros((3,3)))
wyniki["Cy_cor_iter3"].append(numpy.zeros((3,3)))


wyniki["Cx"].append(0.2) 
wyniki["Cy"].append(0.2)
wyniki["psi_0"].append(numpy.array([[0.,0.,0.,],[0.,1.,0.],[0.,0.,0.,]]))
wyniki["psi_iter1"].append(numpy.array([[0.,0.,0.,],[0.,.6,.2],[0.,.2,0.,]]))
wyniki["psi_iter2"].append(numpy.array([[0.,0.,0.,],[0.,.64,.18],[0.,.18,0.,]]))
wyniki["psi_iter3"].append(numpy.array([[0.,0.,0.,],[0.,.657,.171],[0.,.171,0.,]]))
wyniki["Ax_iter2"].append(numpy.array([[0.,1.,1.],[0.,-.5,-1.],[0.,-1.,0.]]))
wyniki["Ay_iter2"].append(wyniki["Ax_iter2"][-1].transpose())
wyniki["Ax_iter3"].append(numpy.array([[0.,1.,1.,],[0.,-.5610,-1.],[0.,-1.,0.,]]))
wyniki["Ay_iter3"].append(wyniki["Ax_iter3"][-1].transpose())
wyniki["Bx_iter2"].append(0.5*numpy.array([[0.5,1.,-1.,],[0.6,1.,-1.],[1.,0.,-1.,]]))
wyniki["By_iter2"].append(wyniki["Bx_iter2"][-1].transpose())
wyniki["Bx_iter3"].append(0.5*numpy.array([[0.5610,1.,-1.,],[0.64,1.,-1.],[1.,0.,-1.,]]))
wyniki["By_iter3"].append(wyniki["Bx_iter3"][-1].transpose())
wyniki["Cx_cor_iter2"].append(numpy.array([[-.02,-.16,.02,],[-.01,.14,.18],[-.012,-.1,-.14,]]))
wyniki["Cy_cor_iter2"].append(wyniki["Cx_cor_iter2"][-1].transpose())
wyniki["Cx_cor_iter3"].append(numpy.array([[4.2e-4,-.1344,-3.8e-4,],[-.3506,.1201,.14103],[1.728e-4,-4.84e-2,-.11375]]))
wyniki["Cy_cor_iter3"].append(wyniki["Cx_cor_iter3"][-1].transpose())


wyniki["Cx"].append(0.5) 
wyniki["Cy"].append(0.5)
wyniki["psi_0"].append(numpy.array([[0.,0.,0.,],[0.,1.,0.],[0.,0.,0.,]]))
wyniki["psi_iter1"].append(numpy.array([[0.,0.,0.,],[0.,0.,.5],[0.,.5,0.,]]))
wyniki["psi_iter2"].append(wyniki["psi_iter1"][-1])
wyniki["psi_iter3"].append(wyniki["psi_iter1"][-1])
wyniki["Ax_iter2"].append(numpy.array([[0.,0.,1.],[0.,1.,-1.],[0.,-1.,0.]]))
wyniki["Ay_iter2"].append(wyniki["Ax_iter2"][-1].transpose())
wyniki["Ax_iter3"].append(wyniki["Ax_iter2"][-1])
wyniki["Ay_iter3"].append(wyniki["Ax_iter3"][-1].transpose())
wyniki["Bx_iter2"].append(0.5*numpy.array([[-1.,1.,0.,],[0.,1.,-1.],[1.,0.,-1.,]]))
wyniki["By_iter2"].append(wyniki["Bx_iter2"][-1].transpose())
wyniki["Bx_iter3"].append(wyniki["Bx_iter2"][-1])
wyniki["By_iter3"].append(wyniki["Bx_iter3"][-1].transpose())
wyniki["Cx_cor_iter2"].append(numpy.array([[-.125,-.25,.125,],[.125,-.125,.25],[0.,.125,-.125,]]))
wyniki["Cy_cor_iter2"].append(wyniki["Cx_cor_iter2"][-1].transpose())
wyniki["Cx_cor_iter3"].append(numpy.array([[0.005859,-.1875,-.001953,],[-.005859,.001953,.1875],[0.,.10742,-.10742]]))
wyniki["Cy_cor_iter3"].append(wyniki["Cx_cor_iter3"][-1].transpose())


def wypisywanie(wybor_lista):
    for i in wybor_lista:
        print "numer w liscie", i
        for k in wyniki.keys():
            print k, '\n', wyniki[k][i]

wypisywanie([0,1,2,3])
