
__version__ = "0.1"

import itertools
from warnings import warn

from colomoto.minibn import BooleanNetwork

from colomoto.types import *
from algorecell_types import *

from .iface import CabeanIface

from .debug import enable_debug, disable_debug

def alias(a):
    return "a{}".format(a)

class _CabeanReprogramming(object):
    def __init__(self, bn):
        if not isinstance(bn, BooleanNetwork):
            raise TypeError("argument must be a minibn.BooleanNetwork object")
        self.bn = bn

    def register_aliases(self, strategies, attractor_ids):
        for a in attractor_ids:
            strategies.register_alias(alias(a), self.attractors[a])

    def strategy_step(self, a, m, next_step=None):
        orig = self.attractors[a]
        #m = assignments_from_flips(orig, flips)
        p = InstantaneousPerturbation(m)
        t = FromSteadyState if orig.is_single_state else FromOneInLimitCycle
        return t(alias(a), p, *((next_step,) if next_step is not None else ()))


def matching_attractors(attractors, pstate):
    return [i for i,a in attractors.items() if a.match_partial_state(pstate)]

def assignments_from_flips(orig, nodes):
    orig_state = orig.project(nodes)
    return dict([(k,((1-v) if v in [0,1] else "-{}".format(k))) \
            for (k,v) in orig_state.items()])

class _CabeanAttractorReprogramming(_CabeanReprogramming):
    def __init__(self, bn, inputs=None):
        """
        TODO
        """
        assert not inputs or set(bn.inputs()).issuperset(inputs.keys()),\
                "specified inputs are not input nodes of the Boolean network"
        self.iface = CabeanIface(bn, init=inputs)
        self.attractors = self.iface.attractors()

    def check_attractors_integrity(self, result, *indexes):
        if not indexes:
            indexes = self.attractors.keys()
        for i in indexes:
            if self.attractors[i] != result.attractors[i]:
                warn("CABEAN: unstable indexes of attractors... trying again...")
                self.attractors = result.attractors
                return False
        return True

class OneStepReprogramming(_CabeanAttractorReprogramming):
    """
    TODO
    """
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
                        "-control", "OI", "-sin", str(a+1), "-tin", str(b+1),
                        *args)
                if not self.check_attractors_integrity(result, a, b):
                    return self.attractor_to_attractor(orig, dest)
                controls = result.parse_OI()
                for sol in controls.get((a,b),[]):
                    s = self.strategy_step(a, sol)
                    strategies.add(s, result=alias(b))
        return strategies

class AttractorSequentialReprogramming(_CabeanAttractorReprogramming):
    """
    TODO
    """
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
        used_attractors = set(aorigs).union(adests)
        for a in aorigs:
            for b in adests:
                result = self.iface.execute("-compositional", "2",
                        "-control", "ASI", "-sin", str(a+1), "-tin", str(b+1),
                        *args)
                if not self.check_attractors_integrity(result, a, b):
                    return self.attractor_to_attractor(orig, dest)
                controls = result.parse_ASI()
                for sol in controls.get((a,b),[]):
                    s = None
                    for (c, m) in reversed(sol):
                        s = self.strategy_step(c, m, s)
                        used_attractors.add(c)
                    strategies.add(s, result=alias(b))
        self.register_aliases(strategies, used_attractors)
        return strategies

class SequentialReprogramming(_CabeanReprogramming):
    def attractor_to_attractor(self, orig, dest, maxsteps=5, limit=200):
        if limit == 1:
            l = "1"
        elif limit <= 200:
            l = "0"
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

def attractors(bn, *spec, **kwspec):
    """
    TODO
    """
    init = PartialState(*spec, **kwspec)
    assert set(bn.inputs()).issuperset(init.keys()),\
            "specified inputs are not input nodes of the Boolean network"
    iface = CabeanIface(bn, init=init)
    return list(iface.attractors().values())

__all__ = ["attractors",
        "OneStepReprogramming",
        "SequentialReprogramming",
        "AttractorSequentialReprogramming"]
