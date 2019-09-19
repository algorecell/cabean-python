
__version__ = "0.1"

import itertools

from colomoto.minibn import BooleanNetwork

from colomoto.reprogramming import *
from colomoto.types import *

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
            strategies.register_alias(alias(a), self.result.attractors[a])

    def strategy_step(self, a, flips, next_step=None):
        orig = self.result.attractors[a]
        m = assignments_from_flips(orig, flips)
        p = InstantaneousPerturbation(m)
        t = FromSteadyState if orig.is_single_state else FromOneInLimitCycle
        return t(alias(a), p, *((next_step,) if next_step is not None else ()))


def matching_attractors(attractors, pstate):
    return [i for i,a in attractors.items() if a.match_partial_state(pstate)]

def assignments_from_flips(orig, nodes):
    orig_state = orig.project(nodes)
    return dict([(k,((1-v) if v in [0,1] else "-{}".format(k))) \
            for (k,v) in orig_state.items()])



class OneStepReprogramming(_CabeanReprogramming):
    """
    TODO
    """
    def __init__(self, bn, inputs=None):
        """
        TODO
        """
        assert not inputs or set(bn.inputs()).issuperset(inputs.keys(),
            "specified inputs are not input nodes of the Boolean network"
        self.iface = CabeanIface(bn, init=inputs)
        self.result = self.iface.execute("-compositional", "2", "-control", "1")
        self.controls = self.result.parse_control1()

    def attractor_to_attractor(self, orig, dest):
        """
        TODO
        """
        aorigs = matching_attractors(self.result.attractors, orig)
        adests = matching_attractors(self.result.attractors, dest)
        strategies = ReprogrammingStrategies()
        self.register_aliases(strategies, set(aorigs).union(adests))
        for a in aorigs:
            for b in adests:
                for sol in self.controls.get((a,b),[]):
                    s = self.strategy_step(a, sol)
                    strategies.add(s, result=alias(b))
        return strategies

class AttractorSequentialReprogramming(_CabeanReprogramming):
    """
    TODO
    """
    def __init__(self, bn, inputs=None):
        """
        TODO
        """
        assert not inputs or set(bn.inputs()).issuperset(inputs.keys(),
            "specified inputs are not input nodes of the Boolean network"
        self.iface = CabeanIface(bn, init=inputs)
        self.result = self.iface.execute("-compositional", "2", "-control", "3")
        self.controls = self.result.parse_control3()

    def attractor_to_attractor(self, orig, dest):
        """
        TODO
        """
        aorigs = matching_attractors(self.result.attractors, orig)
        adests = matching_attractors(self.result.attractors, dest)
        strategies = ReprogrammingStrategies()
        used_attractors = set(aorigs).union(adests)
        for a in aorigs:
            for b in adests:
                for sol in self.controls.get((a,b),[]):
                    s = None
                    for (c, flips) in reversed(sol):
                        s = self.strategy_step(c, flips, s)
                        used_attractors.add(c)
                    strategies.add(s, result=alias(b))
        self.register_aliases(strategies, used_attractors)
        return strategies

class SequentialReprogramming(_CabeanReprogramming):
    def attractor_to_attractor(self, orig, dest, maxsteps=5, limit=200):
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
        result = iface.execute("-control", "2", "-path", l)
        controls = result.parse_control2()
        for steps in itertools.product(*controls):
            s = None
            for step in reversed(steps):
                m = assignments_from_flips(State(step["from"]), step["flip"])
                p = InstantaneousPerturbation(m)
                sa = state_alias(step["from"])
                s = FromState(sa, p, *((s,) if s is not None else ()))
            strategies.add(s)
        return strategies

def attractors(bn, *spec, **kwspec):
    """
    TODO
    """
    init = PartialState(*spec, **kwspec)
    assert set(bn.inputs()).issuperset(init.keys(),
            "specified inputs are not input nodes of the Boolean network"
    iface = CabeanIface(bn, init=init)
    result = iface.execute("-compositional", "2")
    return list(result.parse_attractors().values())

__all__ = ["attractors",
        "OneStepReprogramming",
        "SequentialReprogramming",
        "AttractorSequentialReprogramming"]
