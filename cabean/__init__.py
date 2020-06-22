"""
TODO
"""

import itertools
from warnings import warn
import sys

from colomoto.minibn import BooleanNetwork

from colomoto.types import *
from algorecell_types import *

from .iface import CabeanIface

from .debug import *

def load(bn, *spec, **kwspec):
    """
    TODO

    :rtype: :py:class:`.CabeanInstance`
    """
    return CabeanInstance(bn, *spec, **kwspec)

class CabeanInstance(object):
    """
    TODO
    """
    def __init__(self, bn, *spec, **kwspec):
        bn = BooleanNetwork.auto_cast(bn)
        init = PartialState(*spec, **kwspec)
        assert set(bn.inputs()).issuperset(init.keys()),\
                "specified inputs are not input nodes of the Boolean network"
        self.iface = CabeanIface(bn, init=init)
        self.attractors = self.iface.attractors()

def _cabean_instance(model, *spec, **kwspec):
    if not isinstance(model, CabeanInstance):
        return CabeanInstance(model, *spec, **kwspec)
    if spec or kwspec:
        raise TypeError("wrong arguments")
    return model

def alias(a):
    return "a{}".format(a)

_PTYPE = {
    "I": InstantaneousPerturbation,
    "T": TemporaryPerturbation,
    "P": PermanentPerturbation,
}

def assignments_from_flips(orig, nodes):
    orig_state = orig.project(nodes)
    return dict([(k,((1-v) if v in [0,1] else "-{}".format(k))) \
            for (k,v) in orig_state.items()])

class _CabeanReprogramming(object):

    def register_aliases(self, strategies, attractor_ids):
        for a in attractor_ids:
            strategies.register_alias(alias(a), self.attractors[a])

    def strategy_step(self, a, m, next_step=None):
        orig = self.attractors[a]
        p = _PTYPE[self.method[-1]](m)
        t = FromSteadyState if orig.is_single_state else FromOneInLimitCycle
        return t(alias(a), p, *((next_step,) if next_step is not None else ()))


def matching_attractors(attractors, pstate):
    return [i for i,a in attractors.items() if a.match_partial_state(pstate)]

class _CabeanAttractorReprogramming(_CabeanReprogramming):
    def __init__(self, bn, inputs=None):
        self.ci = _cabean_instance(bn, inputs) if inputs else _cabean_instance(bn)
        self.iface = self.ci.iface
        self.attractors = self.ci.attractors

    def check_attractors_integrity(self, result, *indexes):
        if debug_enabled():
            if not indexes:
                indexes = self.attractors.keys()
            for i in indexes:
                if self.attractors[i] != result.attractors[i]:
                    warn("CABEAN: unstable indexes of attractors... trying again...")
                    self.attractors = result.attractors
                    return False
            return True
        return True

class _OneStep(_CabeanAttractorReprogramming):
    def attractor_to_attractor(self, orig, dest, exclude=None):
        """
        TODO
        """
        args = []
        if exclude:
            args += ["-rmPert",
                    self.iface.make_exclude_perturbations(exclude)]
        aorigs = matching_attractors(self.attractors, orig)
        adests = matching_attractors(self.attractors, dest)
        strategies = ReprogrammingStrategies()
        self.register_aliases(strategies, set(aorigs).union(adests))
        for a in aorigs:
            for b in adests:
                result = self.iface.execute("-compositional", "2",
                        "-control", self.method, "-sin", str(a+1), "-tin", str(b+1),
                        *args)
                if not self.check_attractors_integrity(result, a, b):
                    return self.attractor_to_attractor(orig, dest)
                controls = getattr(result, f"parse_{self.method}")()
                for sol in controls.get((a,b),[]):
                    s = self.strategy_step(a, sol)
                    strategies.add(s, result=alias(b))
        return strategies

class OneStep_Instantaneous(_OneStep):
    """
    TODO
    """
    method = "OI"
class OneStep_Temporary(_OneStep):
    """
    TODO
    """
    method = "OT"
class OneStep_Permanent(_OneStep):
    """
    TODO
    """
    method = "OP"

class _AttractorSequential(_CabeanAttractorReprogramming):
    """
    TODO
    """
    def attractor_to_attractor(self, orig, dest, exclude=None, maxpert=None):
        """
        TODO
        """
        args = []
        if exclude:
            args += ["-rmPert",
                    self.iface.make_exclude_perturbations(exclude)]
        if maxpert:
            args += ["-maxpert", str(maxpert)]
        aorigs = matching_attractors(self.attractors, orig)
        adests = matching_attractors(self.attractors, dest)
        strategies = ReprogrammingStrategies()
        used_attractors = set(aorigs).union(adests)
        for a in aorigs:
            for b in adests:
                result = self.iface.execute("-compositional", "2",
                        "-control", self.method, "-sin", str(a+1), "-tin", str(b+1),
                        *args)
                if not self.check_attractors_integrity(result, a, b):
                    return self.attractor_to_attractor(orig, dest)
                controls = getattr(result, f"parse_{self.method}")()
                for sol in controls.get((a,b),[]):
                    s = None
                    for (c, m) in reversed(sol):
                        s = self.strategy_step(c, m, s)
                        used_attractors.add(c)
                    strategies.add(s, result=alias(b))
        self.register_aliases(strategies, used_attractors)
        return strategies

class AttractorSequential_Instantaneous(_AttractorSequential):
    """
    TODO
    """
    method = "ASI"
class AttractorSequential_Temporary(_AttractorSequential):
    """
    TODO
    """
    method = "AST"
class AttractorSequential_Permanent(_AttractorSequential):
    """
    TODO
    """
    method = "ASP"

class Sequential_Instantaneous(_CabeanReprogramming):
    """
    TODO
    """
    method ="GSI"
    def __init__(self, bn):
        if isinstance(bn, CabeanInstance):
            bn = bn.bn
        else:
            bn = BooleanNetwork.auto_cast(bn)
        self.bn = bn

    def attractor_to_attractor(self, orig, dest, maxsteps=5, limit=1):
        """
        TODO
        """
        if limit == 1:
            l = "1"
        else:
            l = "2"

        strategies = ReprogrammingStrategies()

        state2alias = {}
        def state_alias(s):
            frozen = tuple(sorted(s.items()))
            sa = state2alias.get(frozen, None)
            if sa is None:
                sa = "s{}".format(len(state2alias))
                state2alias[frozen] = sa
                strategies.register_alias(sa, s)
            return sa

        iface = CabeanIface(self.bn, pc=maxsteps, init=orig, red=dest)
        result = iface.execute("-control", "GSI", "-path", l)
        controls = result.parse_GSI()
        if not controls:
            return strategies
        i = 0
        for steps in itertools.product(*controls):
            s = None
            for step in reversed(steps):
                m = assignments_from_flips(State(step["from"]), step["flip"])
                p = InstantaneousPerturbation(m)
                sa = state_alias(step["from"])
                s = FromState(sa, p, *((s,) if s is not None else ()))
            strategies.add(s)
            i += 1
            if i == limit:
                break
        return strategies

def attractors(model, *spec, **kwspec):
    """
    TODO
    """
    ci = _cabean_instance(model, *spec, **kwspec)
    return list(ci.attractors.values())

