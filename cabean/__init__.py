"""
This Python module is an interface to the software tool CABEAN
(https://satoss.uni.lu/software/CABEAN/) for the computation and source-target
control of attractors of asynchronous Boolean networks, using symbolic encoding
of state transition graph.

The control predictions can be processed using the `algorecell_types
<https://github.com/algorecell/algorecell_types>`_ library, which eases the
display and comparison with other control methods.

Installation instructions at https://github.com/algorecell/cabean-python.

Examples can be found at:
    - https://nbviewer.jupyter.org/github/algorecell/cabean-python/tree/main/examples/

Quick usage:

>>> import cabean

Model loading:

>>> cb = cabean.load("model.bnet") # in BoolNet format
# alternatively, load with biolqm in any format
>>> import biolqm
>>> lm = biolqm.load("model.zginml") # or any format support by bioLQM
>>> cb = cabean.load(lm)

Attractors :

>>> from colomoto_jupyter import tabulate
>>> tabulate(cabean.attractors(cb))

Reprogramming predictions:

>>> source = {"A": 0, "B": 1} # any attractor matching
>>> target = {"A": 1, "B", 0} # any attractor matching
# one-step reprogramming with instantaneeous perturbations:
>>> cb_inst = cabean.OneStep_Instantaneous(cb)
>>> rs = cb_inst.attractor_to_attractor(source, target)
>>> rs.as_table()
# attractor-sequential reprogramming with temporary perturbations:
>>> cb_seq_temp = cabean.AttractorSequential_Temporary(cb)
>>> rs = cb_seq_temp.attractor_to_attractor(source, target)
>>> rs.as_graph()

See ``help(rs)`` for other display methods

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
    Returns :py:class:`.CabeanInstance` `(bn, *spec, **kwspec)`

    Example :

    >>> cb = cabean.load(bn, {"I1": 1, "I2": 0})
    """
    return CabeanInstance(bn, *spec, **kwspec)

class CabeanInstance(object):
    """
    CABEAN Boolean network model, storing the list of its attractors
    """
    def __init__(self, bn, *spec, **kwspec):
        """
        :param bn: Boolean network in any format supported by
            ``colomoto.minibn.BooleanNetwork``, which include filename in BoolNet
            format, and ``biolqm`` or ``ginsim`` objects.

        Other arguments can used to specify a fixed value for input nodes.

        Example :

        >>> cb = cabean.CabeanInstance(bn, {"I1": 1, "I2": 0})
        """
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
    """
    One-step reprogramming strategies consist of a set of perturbations which,
    when applied in the initial state are sufficient to ensure the reachability
    of the target attractor.
    """
    def attractor_to_attractor(self, orig, dest, exclude=None):
        """
        Compute one-step reprogramming strategies for enforcing the reachability
        of an attractor of the model matching with `dest` from an attractor matching with
        `orig`.
        The type of perturbations depends on the class from which this method is
        called:

        * :py:class:`.OneStep_Instantaneous`
        * :py:class:`.OneStep_Temporary`
        * :py:class:`.OneStep_Permanent`

        :keyword list(str) exclude: list of nodes to exclude from perturbations.

        :rtype: `algorecell_types.ReprogrammingStrategies <https://algorecell-types.readthedocs.io/#algorecell_types.ReprogrammingStrategies>`_
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
    One-step reprogramming with instantaneous perturbations
    """
    method = "OI"
class OneStep_Temporary(_OneStep):
    """
    One-step reprogramming with temporary perturbations, i.e., to be released
    once in an attractor
    """
    method = "OT"
class OneStep_Permanent(_OneStep):
    """
    One-step reprogramming with permanent perturbations. Note that CABEAN only
    considers attractors of the "wild-type" Boolean networks: permanent
    reprogramming strategies correspond to the cases whenever temporary
    perturbations can be sustained forever.
    """
    method = "OP"

class _AttractorSequential(_CabeanAttractorReprogramming):
    """
    Attractor-sequential reprogramming strategies consider reprogramming in
    several steps, where perturbations are performed once in an attractor, and
    may go through several intermediate attractors before reaching the target
    one.
    """
    def attractor_to_attractor(self, orig, dest, exclude=None, maxpert=None):
        """
        Compute attractor-sequential reprogramming strategies for enforcing the reachability
        of an attractor of the model matching with `dest` from an attractor matching with
        `orig`.
        The type of perturbations depends on the class from which this method is
        called:

        * :py:class:`.AttractorSequential_Instantaneous`
        * :py:class:`.AttractorSequential_Temporary`
        * :py:class:`.AttractorSequential_Permanent`

        :keyword list(str) exclude: list of nodes to exclude from perturbations.
        :keyword int maxpert: maximum number of steps

        :rtype: `algorecell_types.ReprogrammingStrategies <https://algorecell-types.readthedocs.io/#algorecell_types.ReprogrammingStrategies>`_
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
    Attractor-sequential reprogramming with instantaneous perturbations.
    """
    method = "ASI"
class AttractorSequential_Temporary(_AttractorSequential):
    """
    Attractor-sequential reprogramming with temporary perturbations, i.e., to be released
    once in an attractor
    """
    method = "AST"
class AttractorSequential_Permanent(_AttractorSequential):
    """
    Attractor-sequential reprogramming with permanent perturbations. Note that CABEAN only
    considers attractors of the "wild-type" Boolean networks: permanent
    reprogramming strategies correspond to the cases whenever temporary
    perturbations can be sustained forever.
    """
    method = "ASP"

class Sequential_Instantaneous(_CabeanReprogramming):
    """
    Sequential reprogramming strategies consider reprogramming in several steps,
    where perturbations are performed in specific states, not necessarily in an
    attractor.
    The current implementation only considers instantaneous perturbations.
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
        Compute sequential reprogramming strategies for enforcing the reachability
        of an attractor of the model matching with `dest` from states matching
        `orig`, using instantenous perturbations.

        :keyword list(str) exclude: list of nodes to exclude from perturbations.
        :keyword int maxsteps: maximum number of steps
        :keyword int limit: maximum number of solutions

        :rtype: `algorecell_types.ReprogrammingStrategies <https://algorecell-types.readthedocs.io/#algorecell_types.ReprogrammingStrategies>`_
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
    Returns the list of attractors of `model`.
    An attractor is either represented as a
    ``colomoto.types.TrapSpaceAttractor``, i.e., a dictionnary mapping node
    names to ``0``, ``1``, or ``*`` (cyclic); or as a
    ``colomoto.types.TrapSpacesAttractor``, i.e., a list of the former type.

    Arguments are the same as :py:class:`.CabeanInstance`.
    """
    ci = _cabean_instance(model, *spec, **kwspec)
    return list(ci.attractors.values())

