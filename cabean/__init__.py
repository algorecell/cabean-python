
__version__ = "0.1"

from colomoto.minibn import BooleanNetwork

from colomoto.reprogramming import *
from colomoto.types import *

from .iface import CabeanIface

from .debug import enable_debug, disable_debug

class _CabeanReprogramming(object):
    def __init__(self, bn):
        if not isinstance(bn, BooleanNetwork):
            raise TypeError("argument must be a minibn.BooleanNetwork object")
        self.bn = bn

def matching_attractors(attractors, pstate):
    return [i for i,a in attractors.items() if a.match_partial_state(pstate)]

class OneStepReprogramming(_CabeanReprogramming):
    def __init__(self, bn, init=None):
        self.iface = CabeanIface(bn, init=init)
        self.result = self.iface.execute("-compositional", "2", "-control", "1")
        self.controls = self.result.parse_control1()

    def attractor_to_attractor(self, orig, dest):
        aorigs = matching_attractors(self.result.attractors, orig)
        adests = matching_attractors(self.result.attractors, dest)
        strategies = [] # TODO object to store state aliases
        for a in aorigs:
            for b in adests:
                sols = self.controls.get((a,b),[])
                strategies.append(sols)
        return strategies

class AttractorSequentialReprogramming(_CabeanReprogramming):
    def __init__(self, bn, init=None):
        self.iface = CabeanIface(bn, init=init)
        self.result = self.iface.execute("-compositional", "2", "-control", "3")
        self.controls = self.result.parse_control3()

    def attractor_to_attractor(self, orig, dest):
        aorigs = matching_attractors(self.result.attractors, orig)
        adests = matching_attractors(self.result.attractors, dest)
        strategies = [] # TODO object to store state aliases
        for a in aorigs:
            for b in adests:
                sols = self.controls.get((a,b),[])
                strategies.append(sols)
        return strategies

class SequentialReprogramming(_CabeanReprogramming):
    def attractor_to_attractor(self, orig, dest, maxsteps=5, limit=200):
        if limit == 0 or limit > 200:
            l = "2"
        elif limit == 1:
            l = "1"
        else:
            l = "0"
        controls = []
        if isinstance(orig, (TrapSpaceAttractor, TrapSpacesAttractor)) \
                or isinstance(orig, dict) and len(orig) == len(self.bn):
            aorigs = [orig]
        else:
            aorigs = matching_attractors(attractors(self.bn), orig)
        for aorig in aorigs:
            iface = CabeanIface(self.bn, pc=maxsteps, init=aorig, red=dest)
            result = iface.execute("-control", "2", "-path", l)
            controls.append(result.parse_control2())
        return controls

def attractors(bn, *spec, **kwspec):
    init = PartialState(*spec, **kwspec)
    iface = CabeanIface(bn, init=init)
    result = iface.execute("-compositional", "2")
    return list(result.parse_attractors().values())

__all__ = ["attractors",
        "OneStepReprogramming",
        "SequentialReprogramming",
        "AttractorSequentialReprogramming"]
