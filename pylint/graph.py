# Copyright (c) 2015-2018, 2020 Claudiu Popa <pcmanticore@gmail.com>
# Copyright (c) 2015 Florian Bruhin <me@the-compiler.org>
# Copyright (c) 2016 Ashley Whetter <ashley@awhetter.co.uk>
# Copyright (c) 2018 ssolanki <sushobhitsolanki@gmail.com>
# Copyright (c) 2019, 2021 Pierre Sassoulas <pierre.sassoulas@gmail.com>
# Copyright (c) 2019 Nick Drozd <nicholasdrozd@gmail.com>
# Copyright (c) 2020 hippo91 <guillaume.peillex@gmail.com>
# Copyright (c) 2020 Damien Baty <damien.baty@polyconseil.fr>
# Copyright (c) 2020 谭九鼎 <109224573@qq.com>
# Copyright (c) 2020 Benjamin Graham <benwilliamgraham@gmail.com>
# Copyright (c) 2021 Daniël van Noord <13665637+DanielNoord@users.noreply.github.com>
# Copyright (c) 2021 Marc Mueller <30130371+cdce8p@users.noreply.github.com>
# Copyright (c) 2021 Andreas Finkler <andi.finkler@gmail.com>
# Copyright (c) 2021 Andrew Howe <howeaj@users.noreply.github.com>

# Licensed under the GPL: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
# For details: https://github.com/PyCQA/pylint/blob/main/LICENSE

"""Graph manipulation utilities.

(dot generation adapted from pypy/translator/tool/make_dot.py)
"""
import codecs
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Optional


def target_info_from_filename(filename):
    """Transforms /some/path/foo.png into ('/some/path', 'foo.png', 'png')."""
    basename = os.path.basename(filename)
    storedir = os.path.dirname(os.path.abspath(filename))
    target = os.path.splitext(filename)[-1][1:]
    return storedir, basename, target


class DotBackend:
    """Dot File backend."""

    def __init__(
        self,
        graphname,
        rankdir=None,
        size=None,
        ratio=None,
        charset="utf-8",
        renderer="dot",
        additional_param=None,
    ):
        if additional_param is None:
            additional_param = {}
        self.graphname = graphname
        self.renderer = renderer
        self.lines = []
        self._source = None
        self.emit(f"digraph {normalize_node_id(graphname)} {{")
        if rankdir:
            self.emit(f"rankdir={rankdir}")
        if ratio:
            self.emit(f"ratio={ratio}")
        if size:
            self.emit(f'size="{size}"')
        if charset:
            assert charset.lower() in {
                "utf-8",
                "iso-8859-1",
                "latin1",
            }, f"unsupported charset {charset}"
            self.emit(f'charset="{charset}"')
        for param in additional_param.items():
            self.emit("=".join(param))

    def get_source(self):
        """Returns self._source."""
        if self._source is None:
            self.emit("}\n")
            self._source = "\n".join(self.lines)
            del self.lines
        return self._source

    source = property(get_source)

    def generate(
        self, outputfile: Optional[str] = None, mapfile: Optional[str] = None
    ) -> str:
        """Generates a graph file.

        :param str outputfile: filename and path [defaults to graphname.png]
        :param str mapfile: filename and path

        :rtype: str
        :return: a path to the generated file
        :raises RuntimeError: if the executable for rendering was not found
        """
        graphviz_extensions = ("dot", "gv")
        name = self.graphname
        if outputfile is None:
            target = "png"
            pdot, dot_sourcepath = tempfile.mkstemp(".gv", name)
            ppng, outputfile = tempfile.mkstemp(".png", name)
            os.close(pdot)
            os.close(ppng)
        else:
            _, _, target = target_info_from_filename(outputfile)
            if not target:
                target = "png"
                outputfile = outputfile + "." + target
            if target not in graphviz_extensions:
                pdot, dot_sourcepath = tempfile.mkstemp(".gv", name)
                os.close(pdot)
            else:
                dot_sourcepath = outputfile
        with codecs.open(dot_sourcepath, "w", encoding="utf8") as file:
            file.write(self.source)
        if target not in graphviz_extensions:
            if shutil.which(self.renderer) is None:
                raise RuntimeError(
                    f"Cannot generate `{outputfile}` because '{self.renderer}' "
                    "executable not found. Install graphviz, or specify a `.gv` "
                    "outputfile to produce the DOT source code."
                )
            use_shell = sys.platform == "win32"
            if mapfile:
                subprocess.call(
                    [
                        self.renderer,
                        "-Tcmapx",
                        "-o",
                        mapfile,
                        "-T",
                        target,
                        dot_sourcepath,
                        "-o",
                        outputfile,
                    ],
                    shell=use_shell,
                )
            else:
                subprocess.call(
                    [self.renderer, "-T", target, dot_sourcepath, "-o", outputfile],
                    shell=use_shell,
                )
            os.unlink(dot_sourcepath)
        return outputfile

    def emit(self, line):
        """Adds <line> to final output."""
        self.lines.append(line)

    def emit_edge(self, name1, name2, **props):
        """Emit an edge from <name1> to <name2>.
        edge properties: see https://www.graphviz.org/doc/info/attrs.html
        """
        attrs = [f'{prop}="{value}"' for prop, value in props.items()]
        n_from, n_to = normalize_node_id(name1), normalize_node_id(name2)
        self.emit(f"{n_from} -> {n_to} [{', '.join(sorted(attrs))}];")

    def emit_node(self, name, **props):
        """Emit a node with given properties.
        node properties: see https://www.graphviz.org/doc/info/attrs.html
        """
        attrs = [f'{prop}="{value}"' for prop, value in props.items()]
        self.emit(f"{normalize_node_id(name)} [{', '.join(sorted(attrs))}];")


def normalize_node_id(nid):
    """Returns a suitable DOT node id for `nid`."""
    return f'"{nid}"'


def get_cycles(graph_dict, vertices=None):
    """Given a dictionary representing an ordered graph (i.e. key are vertices
    and values is a list of destination vertices representing edges), return a
    list of detected cycles
    """
    if not graph_dict:
        return ()
    result = []
    if vertices is None:
        vertices = graph_dict.keys()
    for vertice in vertices:
        _get_cycles(graph_dict, [], set(), result, vertice)
    return result

def _get_cycles(graph_dict, path, visited, result, vertice):
    """Recursive function doing the real work for get_cycles."""
    if vertice in path:
        cycle = [str(vertice)]
        for node in path[::-1]:
            if node == vertice:
                break
            cycle.insert(0, str(node))
        # make a canonical representation
        start_from = min(cycle)
        index = cycle.index(start_from)
        cycle = cycle[index:] + cycle[0:index]
        # append it to result if not already in
        if cycle not in result:
            result.append(cycle)
        return
    path.append(vertice)
    try:
        for node in graph_dict[vertice]:
            # don't check already visited nodes again
            if node not in visited:
                _get_cycles(graph_dict, path, visited, result, node)
                visited.add(node)
    except KeyError:
        pass
    path.pop()

# symbol_dict: Dict[Tuple[str, str], str]
# Where the value is ">=", ">", "<=", "<"
def handle_cycles(graph_dict, symbol_dict, cycles):
    for cycle in cycles:
        all_geq = all(symbol_dict[(cycle[i], cycle[i+1])] == ">=" for (i, _) in enumerate(cycle) if i < len(cycle) - 1)
        all_geq = all_geq and symbol_dict[cycle[-1], cycle[0]] == ">="
        if all_geq:
            print("Pylint simplify")
        else:
            print("Pylint bad cycle")

def get_paths(graph_dict, indegree_dict, frequency_dict):
    indegree_zero = sorted([node for node in indegree_dict if indegree_dict[node] == 0], reverse=False, key=str)
    paths = set()
    while indegree_zero:
        path = []
        get_path(path, graph_dict, indegree_zero.pop(), indegree_dict, indegree_zero, frequency_dict)
        print(path)
        if not all(isinstance(item, (int, float)) for item in path):

            # Backtrack to remove redundant links
            for i, item in enumerate(path):
                if isinstance(item, str):
                # if True:
                    for val in path[:i]:
                        graph_dict[val].discard(item)

            path = strip_path(path)
            if len(path) > 1:
                paths.add(tuple(path))
    return paths

def get_path(path, graph_dict, node, indegree_dict, indegree_zero, frequency_dict):
    path.append(node)
    indegree_dict[node] = max(indegree_dict[node] - 1, 0)

    # Sorting makes the results deterministic, also, always visit variable nodes before numerical ones
    adj = sorted([a for a in graph_dict[node] if frequency_dict[(node, a)] != 0], reverse=False, key=str)
    if indegree_dict[node] == 0 and len(adj) >= 2:
        copied = list(path)
        cur = copied.pop()
        while copied and isinstance(copied[-1], (int, float)):
        # while copied:
            cur = copied.pop()
        indegree_zero.append(cur) # This will visit adj[1] when we pop it

    # if len(adj) >= 1 and not indegree_dict[adj[0]] == 0:
    if len(adj) >= 1:
        # graph_dict[node].remove(adj[0])
        frequency_dict[(node, adj[0])] -= 1
        get_path(path, graph_dict, adj[0], indegree_dict, indegree_zero, frequency_dict)

def get_paths_longest(graph_dict, indegree_dict, indegree_zero, frequency_dict):
    paths = []
    while indegree_zero:
        symbols_in_longest_path = {key: 1 if isinstance(key, str) else 0 for key in graph_dict}
        nodes_in_longest_path = {key: 1 for key in graph_dict}
        for root in indegree_zero:
            count_nodes(root, graph_dict, symbols_in_longest_path, nodes_in_longest_path, frequency_dict)
        path = []
        longest_path_item = get_longest_path_item(indegree_zero, symbols_in_longest_path, nodes_in_longest_path)
        indegree_zero.remove(longest_path_item)
        get_path_longest(path, graph_dict, longest_path_item, indegree_dict, indegree_zero, frequency_dict, symbols_in_longest_path, nodes_in_longest_path)
        for i, item in enumerate(path):
            if True:
                for val in path[:i]:
                    # graph_dict[val].discard(item)
                    frequency_dict[(val, item)] = max(frequency_dict[(val, item)] - 1, 0)
        paths.append(path)
        print(path)

    return paths

def get_longest_path_item(items, symbols_in_longest_path, nodes_in_longest_path):
    return sorted(items, key=lambda x: (symbols_in_longest_path[x], nodes_in_longest_path[x]))[0]

def get_path_longest(path, graph_dict, node, indegree_dict, indegree_zero, frequency_dict, symbols_in_longest_path, nodes_in_longest_path):
    path.append(node)
    indegree_dict[node] = max(indegree_dict[node] - 1, 0)
    adj = [a for a in graph_dict[node] if frequency_dict[(node, a)] != 0]
    if len(adj) >= 1:
        next = get_longest_path_item(adj, symbols_in_longest_path, nodes_in_longest_path)
        # frequency_dict[(node, next)] -= 1
        get_path_longest(path, graph_dict, next, indegree_dict, indegree_zero, frequency_dict, symbols_in_longest_path, nodes_in_longest_path)
        if indegree_dict[node] == 0 and (len(adj) >= 2 or frequency_dict[(node, next)] >= 2):
            indegree_zero.add(node)

def count_nodes(node, graph_dict, symbols_in_longest_path, nodes_in_longest_path, frequency_dict):
    for adj in graph_dict[node]:
        if frequency_dict[(node, adj)] == 0:
            continue
        if isinstance(node, str):
            symbols_in_longest_path[adj] = max(symbols_in_longest_path[adj], symbols_in_longest_path[node] + 1)
        nodes_in_longest_path[adj] = max(nodes_in_longest_path[adj], nodes_in_longest_path[node] + 1)
        count_nodes(adj, graph_dict, symbols_in_longest_path, nodes_in_longest_path, frequency_dict)


def strip_path(path):
    low = 0
    high = len(path) - 1
    if path and isinstance(path[0], (int, float)):
        while low < len(path) - 1 and isinstance(path[low], (int, float)) and isinstance(path[low+1], (int, float)):
            low += 1

    if path and isinstance(path[-1], (int, float)):
        while high > 0 and high > low and isinstance(path[high], (int, float)) and isinstance(path[high-1], (int, float)):
            high -= 1
    return path[low:high+1]

import collections
from astroid import nodes
def optimize_boolop(node: nodes.BoolOp):
    if node.op != 'and' or len(node.values) < 2:
        return

    graph_dict = collections.defaultdict(set)
    symbol_dict = {}
    frequency_dict = collections.defaultdict(int)
    indegree_dict = collections.defaultdict(int)
    const_values: list[int] = []

    for statement in node.values:
        ops = list(statement.ops)
        left_statement = statement.left
        while ops:
            if not isinstance(statement, nodes.Compare):
                continue

            left = get_compare_operand_value(left_statement, const_values)
            if left is None:
                continue

            operator, right_statement = ops.pop(0)
            right = get_compare_operand_value(right_statement, const_values)
            if right is None:
                continue

            # Make the graph always point from larger to smaller
            if operator == "<":
                operator = ">"
                left, right = right, left
            elif operator == "<=":
                operator = ">="
                left, right = right, left

            # Update maps
            graph_dict[left].add(right)
            graph_dict[right] # Ensure the node exists in graph
            symbol_dict[(left, right)] = operator
            indegree_dict[left] += 0 # Make sure every node has an entry
            indegree_dict[right] += 1
            frequency_dict[(left, right)] += 1

            left_statement = right_statement

    # Link up constant nodes
    sorted_consts = sorted(const_values)
    while sorted_consts:
        largest = sorted_consts.pop()
        for smaller in set(sorted_consts):
            if smaller < largest:
                symbol_dict[(largest, smaller)] = ">"
                indegree_dict[smaller] += 1
                frequency_dict[(largest, smaller)] += 1
                graph_dict[largest].add(smaller)

                for adj in graph_dict[smaller]:
                    if isinstance(adj, str):
                        graph_dict[largest].discard(adj)

    return graph_dict, symbol_dict, indegree_dict, frequency_dict

def get_compare_operand_value(node: nodes.Compare, const_values: Optional[set[int]]=None):
    value = None
    if isinstance(node, nodes.Name):
        value = node.name
    elif isinstance(node, nodes.Const):
        value = node.value
        const_values.append(value)
    # elif #Walrus
    return value

if __name__ == "__main__":
    import astroid
    node1 = astroid.extract_node("""
    if a < 20 and a < 15 and d < 15 and 15 < e:
        pass
    """)

    node2 = astroid.extract_node("""
    if a < 20 and a < 15 and d < 15 and 15 < e < 20:
        pass
    """)

    node25 = astroid.extract_node("""
    if d < 15 and a < 15 < e < 20:
        pass
    """)


    node3 = astroid.extract_node("""
    if a < 20 and a < 15:
        pass
    """)

    node4 = astroid.extract_node("""
    if a < 20 and a < 15 and b < 20 and b < 15:
        pass
    """)

    node5 = astroid.extract_node("""
    if a < b and a < c < d and c < 20 and b < d and d < f:
        pass
    """)

    node6 = astroid.extract_node("""
    if a < b < c and a < c:
        pass
    """)
    node = node1

    print(node.test.as_string())
    graph_dict, symbol_dict, indegree_dict, frequency_dict = optimize_boolop(node.test)
    print(graph_dict)
    print(symbol_dict)
    print(indegree_dict)
    print(frequency_dict)
    cycles = get_cycles(graph_dict)
    print("cycles", cycles)

    if len(cycles) == 0:
        indegree_zero = {node for node in indegree_dict if indegree_dict[node] == 0}
        paths = get_paths_longest(graph_dict, indegree_dict, indegree_zero, frequency_dict)

        # paths = get_paths(graph_dict, indegree_dict, frequency_dict)
        print("paths", paths)
