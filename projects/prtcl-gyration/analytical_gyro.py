import numpy as np
from py_runko_aux import *

def analytical_gyro_full(time,conf):
    vxref = []
    vyref = []
    xref = []
    yref = []

    # Relativistic gyro frequency (cgs): omega_B = qB/(gamma*m*c)
    # Relativistic gyro frequency (runko): omega_B*dt = qh*Bh*/(gamma*mh*ch) * m-/m+
    lfreq = (conf.qe*conf.binit)/(conf.gamma*abs(conf.qe)*conf.cfl**2)

    # Phase-lag, particles start at a=0 radians
    a = 0

    # Parameter shorthands
    vy = conf.vy
    vy2 = conf.vy2

    cx = conf.NxMesh/2.
    cy = conf.NyMesh/2.

    for t in time:
        vxref.append(np.array([-vy*np.sin(lfreq*t+a),-vy2*np.sin(lfreq*t+a)]))
        vyref.append(np.array([vy*np.cos(lfreq*t+a),vy2*np.cos(lfreq*t+a)]))

        xref.append(np.array([cx+vy/lfreq*np.cos(lfreq*t+a),
                              cx+vy2/lfreq*np.cos(lfreq*t+a)]))

        yref.append(np.array([cy+vy/lfreq*np.sin(lfreq*t+a),
                              cy+vy2/lfreq*np.sin(lfreq*t+a)]))

    xref = np.array(xref)
    yref = np.array(yref)
    vxref = np.array(vxref)
    vyref = np.array(vyref)

    return xref,yref,vxref,vyref


def analytical_gyro_single(time,conf):
    # Relativistic gyro frequency (cgs): omega_B = qB/(gamma*m*c)
    # Relativistic gyro frequency (runko): omega_B*dt = qh*Bh*/(gamma*mh*ch) * m-/m+
    lfreq = (conf.qe*conf.binit)/(conf.gamma*abs(conf.qe)*conf.cfl**2)

    # Phase-lag, particles start at a=0 radians
    a = 0

    # Parameter shorthands
    vy = conf.vy
    vy2 = conf.vy2

    cx = conf.NxMesh/2.
    cy = conf.NyMesh/2.

    t = time

    vxref = np.array([-vy*np.sin(lfreq*t+a),-vy2*np.sin(lfreq*t+a)])
    vyref = np.array([vy*np.cos(lfreq*t+a),vy2*np.cos(lfreq*t+a)])

    xref = np.array([cx+vy/lfreq*np.cos(lfreq*t+a),
                          cx+vy2/lfreq*np.cos(lfreq*t+a)])

    yref = np.array([cy+vy/lfreq*np.sin(lfreq*t+a),
                          cy+vy2/lfreq*np.sin(lfreq*t+a)])


    return xref,yref,vxref,vyref
