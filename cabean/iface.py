import itertools
import os
import subprocess
import tempfile

from colomoto.types import *
from colomoto_jupyter.sessionfiles import new_output_file

from cabean.debug import debug_enabled

cabean_base_options = []

class CabeanProcessError(subprocess.CalledProcessError):
    """
    Exception raised when a Pint command fails.
    """
    def __str__(self):
        stderr = "\n%s" % self.stderr.decode() if self.stderr else self.output.decode()
        return "Command '%s' returned non-zero exit status %d%s" \
            % (" ".join(self.cmd), self.returncode, stderr)

class CabeanResult(object):
    def __init__(self, iface, content):
        self.iface = iface
        self.lines = content.split("\n")
        if debug_enabled():
            print(content)

    @property
    def attractors(self):
        if not hasattr(self, "_CabeanResult__cache_attractors"):
            self.__cache_attractors = self.parse_attractors()
        return self.__cache_attractors

    def parse_state(self, spec):
        spec = spec[0:len(spec):2]
        spec = zip(self.iface.ordered_nodes, spec)
        return PartialState([(x,int(v) if v != "-" else "*") for x,v in spec])

    def parse_attractors(self):
        attractors = {}
        num = None
        for line in self.lines:
            if line.startswith("=") and "=== find attractor #" in line:
                parts = line.split()
                num = int(parts[3][1:])-1
                size = int(parts[5])
            elif num is not None:
                if line.startswith(":"):
                    pass
                elif not line:
                    # TODO: sanity check with size
                    num = None
                else:
                    state = self.parse_state(line.split()[0])
                    state = TrapSpaceAttractor(state)
                    if num not in attractors:
                        attractors[num] = state
                    else:
                        attractors[num] = attractors[num].extend(state)
        return attractors

    def parse_onestep(self, mode):
        controls = {}
        state = 0
        for line in self.lines:
            line = line.strip()
            if f"= ONE-STEP {mode.upper()}" in line:
                state = 1
            elif state == 1 and line.startswith("source -"):
                w = line.split()
                a1, a2 = int(w[2])-1, int(w[5])-1
                controls[(a1,a2)] = []
                state = 2
            elif state == 2 and line.lower().startswith("control set:"):
                line = line[12:]
                p = self.parse_controlset(line)
                controls[(a1,a2)].append(p)
            elif state == 2 and line.startswith("execution time"):
                state = 1
        if debug_enabled():
            print(controls)
        return controls

    def parse_controlset(self, data):
        p = {}
        for c in data.strip().split():
            node, value = c.split("=")
            p[node] = int(value)
        return p

    def parse_attractor_sequential(self, mode):
        controls = {}
        state = 0
        for line in self.lines:
            line = line.strip()
            if line.startswith(f"========= ATTRACTOR-BASED SEQUENTIAL {mode.upper()} "):
                state = 1
            elif state == 1 and line.startswith("source -"):
                w = line.split()
                a1, a2 = int(w[2])-1, int(w[5])-1
                controls[(a1,a2)] = []
            elif state == 1 and line.startswith("Sequence of the attractors"):
                w = line.strip().split()
                seq = [int(aid)-1 for aid in w[4:len(w):2]]
                steps = []
                state = 2
            elif state == 2 and (not line or line.startswith("execution time")):
                for path in itertools.product(*steps):
                    controls[(a1,a2)].append(list(zip(seq[:-1], path)))
                state = 1
            elif state == 2 and line.lower().startswith("step"):
                step = []
                steps.append(step)
            elif state == 2 and line.lower().startswith("control set"):
                line = ":".join(line.split(":")[1:])
                p = self.parse_controlset(line)
                step.append(p)
        if debug_enabled():
            print(controls)
        return controls

    def parse_OI(self):
        return self.parse_onestep("instantaneous")
    def parse_OT(self):
        return self.parse_onestep("temporary")
    def parse_OP(self):
        return self.parse_onestep("permanent")

    def parse_ASI(self):
        return self.parse_attractor_sequential("instantaneous")
    def parse_AST(self):
        return self.parse_attractor_sequential("temporary")
    def parse_ASP(self):
        return self.parse_attractor_sequential("permanent")

    def parse_GSI(self):
        controls = []
        step = 0
        mode = 0
        for line in self.lines:
            line = line.strip()
            if line.startswith("One sequential"):
                mode = 1
            if line.startswith("STEP "):
                step = int(line.split()[1])
                controls.append([])
                assert len(controls) == step
            if line.startswith("path "):
                if mode == 1:
                    controls.append([])
                control = {}
            if line.startswith("from "):
                state = self.parse_state(line.split()[2])
                control["from"] = state
            if line.startswith("driver nodes:"):
                control["flip"] = set(line.split()[2:])
                controls[-1].append(control)
        return controls

    def __str__(self):
        return "\n".join(self.lines)


class CabeanIface(object):
    def __init__(self, bn, init=None, red=None, pc=0):
        constants = bn.constants()
        if constants:
            init = init if init is not None else {}
            bn = bn.copy()
            for n, f in constants.items():
                bn[n] = n
                init[n] = bool(f)
        self.bn = bn
        self.init = init
        self.red = red
        self.pc = pc
        self.ordered_nodes = list(sorted(self.bn.keys()))

    def attractors(self):
        result = self.execute("-compositional", "2")
        return result.attractors

    def make_exclude_perturbations(self, exclude):
        """
        TODO
        """
        excfile = new_output_file(suffix="_rmPert.txt", prefix="cabean")
        x = {"R0": [], "R1": [], "R": []}
        for spec in exclude:
            m = "R"
            n = spec
            if spec.endswith("+"):
                n = spec.strip("+")
                m = "R0"
            elif spec.endswith("-"):
                n = spec.strip("-")
                m = "R1"
            assert n in self.bn, "Unknown node '{}'".format(n)
            x[m].append(n)
        with open(excfile, "w") as fp:
            fp.write("R0: {}\n".format(",".join(x["R0"])))
            fp.write("R1: {}\n".format(",".join(x["R1"])))
            fp.write("R: {}\n".format(",".join(x["R"])))
        return excfile

    def execute(self, *args, isplfile=None):
        args = ["cabean", "-asynbn", "-steadystates", "-newtarjan", "-newpred"] \
                + cabean_base_options + list(args)
        if self.pc:
            args += ["-pc", str(self.pc)]
        if debug_enabled():
            isplfile = "/tmp/cabean.ispl"
        if isplfile:
            tmpfile = isplfile
        else:
            fd, tmpfile = tempfile.mkstemp(suffix=".ispl", prefix="cabean")
            os.close(fd)
        args.append(tmpfile)
        try:
            with open(tmpfile, "w") as fp:
                self.write_ispl(fp)
            output = subprocess.check_output(args, stderr=subprocess.PIPE)
            return CabeanResult(self, output.decode())
        except subprocess.CalledProcessError as e:
            e = CabeanProcessError(e.returncode, e.cmd, e.output, e.stderr)
            raise e from None
        finally:
            if not isplfile:
                os.unlink(tmpfile)

    def write_ispl(self, fp):
        fp.write("Agent M\n\tVars:\n")
        for x in self.ordered_nodes:
            fp.write("\t\t{}: boolean;\n".format(x))
        if self.pc:
            fp.write("\t\tpc: 0..{};\n".format(self.pc))
        fp.write("\tend Vars\n")
        if self.red:
            d = ispl_state(self.red)
            fp.write("\tRedStates:\n\t\t{}\n\tend RedStates\n".format(d))
        fp.write("\tActions = {none};\n")
        fp.write("\tProtocol:\n\t\tOther: {none};\n\tend Protocol\n")
        fp.write("\tEvolution:\n")
        for x in self.ordered_nodes:
            f = self.bn[x]
            if f is self.bn.ba.TRUE:
                fp.write("\t\t{0}=true if {0}=true or {0}=false;\n".format(x))
            elif f is self.bn.ba.FALSE:
                fp.write("\t\t{0}=false if {0}=true or {0}=false;\n".format(x))
            else:
                f = str(f).replace("!", "~")
                fp.write("\t\t{}=true if ({})=true;\n".format(x, f))
                fp.write("\t\t{}=false if ({})=false;\n".format(x, f))
        if self.pc:
            for x in self.ordered_nodes:
                d = {"x": x, "pc": self.pc}
                fp.write("\t\t{x}=true and pc=pc+1 \
                        if pc<{pc} and {x}=false;\n".format(**d))
                fp.write("\t\t{x}=false and pc=pc+1\
                        if pc<{pc} and {x}=true;\n".format(**d))
        fp.write("\tend Evolution\nend Agent\n\n")
        if self.init:
            d = ispl_state(self.init, prefix="M.")
            fp.write("InitStates\n\t{}\nend InitStates\n".format(d))
        else:
            x = self.ordered_nodes[0]
            fp.write("InitStates\n\tM.{0}=true or M.{0}=false;\n\
                        end InitStates\n".format(x))

def ispl_state(state, prefix=""):
    if isinstance(state, list):
        return "\n".join(map(ispl_state, state))
    return "%s;" % " and ".join(["{prefix}{x}={v}".format(x=x,
                    v="true" if v else "false",
                    prefix=prefix) \
                for x,v in state.items() if v != "*"])
