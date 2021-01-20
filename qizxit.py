from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import execute, transpile, Aer

from qiskit.circuit.library import CXGate, TGate, TdgGate, SGate, SdgGate, XGate, ZGate, HGate

from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.dagcircuit import DAGCircuit

from qiskit.tools.visualization import dag_drawer

import pyzx as zx

from functools import reduce

def reduce_t(qc):
    dag = circuit_to_dag(qc)
    fl  = True
    def exchange(gate, pre_gate):
        p    = QuantumRegister(1, "p")
        dag2 = DAGCircuit()
        dag2.add_qreg(p)
        dag2.apply_operation_back(gate,     qargs = p)
        dag2.apply_operation_back(pre_gate, qargs = p)
        dag.remove_op_node(pre_node)
        dag.substitute_node_with_dag(node = node, input_dag = dag2, wires = p)
    def exchange_cx(gate):
        p    = QuantumRegister(2, "p")
        dag2 = DAGCircuit()
        dag2.add_qreg(p)
        dag2.apply_operation_back(gate,     qargs = [p[0]])
        dag2.apply_operation_back(CXGate(), qargs = p)
        dag.remove_op_node(node)
        dag.substitute_node_with_dag(node = pre_node, input_dag = dag2, wires = p)
    def annihilation():
        dag.remove_op_node(node)
        dag.remove_op_node(pre_node)
    def replece(gate):
        p    = QuantumRegister(1, "p")
        dag2 = DAGCircuit()
        dag2.add_qreg(p)
        dag2.apply_operation_back(gate, qargs = p)
        dag.remove_op_node(pre_node)
        dag.substitute_node_with_dag(node = node, input_dag = dag2, wires = p)
        
    while (True):
        fl = True
        # Move T
        for node in dag.op_nodes(op = TGate):
            pre_node = next(dag.quantum_predecessors(node))
            # T + Tdg = annihilation
            if (pre_node.name == 'tdg'):
                annihilation()
                fl = False
            # T + T = S
            if (pre_node.name == 't'):
                replece(SGate())
                fl = False
            # T + S = S + T
            if (pre_node.name == 's'):
                exchange(TGate(), SGate())
                fl = False
            # T + Sdg = Sdg + T
            if (pre_node.name == 'sdg'):
                exchange(TGate(), SdgGate())
                fl = False
            # T + Z = Z + T
            if (pre_node.name == 'z'):
                exchange(TGate(), ZGate())
                fl = False
            # T + CX ctrl = CX ctrl + T
            if (pre_node.name == 'cx' and pre_node.qargs[0] == node.qargs[0]):
                exchange_cx(TGate())
                fl = False
        # Move Tdg
        for node in dag.op_nodes(op = TdgGate):
            pre_node = next(dag.quantum_predecessors(node))
            # Tdg + T = annihilation
            if (pre_node.name == 't'):
                annihilation()
                fl = False
            # Tdg + Tdg = Sdg
            if (pre_node.name == 'tdg'):
                replece(SdgGate())
                fl = False
            # Tdg + S = S + Tdg
            if (pre_node.name == 's'):
                exchange(TdgGate(), Sgate())
                fl = False
            # Tdg + Sdg = Sdg + Tdg
            if (pre_node.name == 'sdg'):
                exchange(TdgGate(), Sdggate())
                fl = False
            # Tdg + Z = Z + Tdg
            if (pre_node.name == 'z'):
                exchange(TdgGate(), Zgate())
                fl = False
            # Tdg + CX ctrl = CX ctrl + Tdg
            if (pre_node.name == 'cx' and pre_node.qargs[0] == node.qargs[0]):
                excahnge_cx(TdgGate())
                fl = False
        # H + H = annihilation
        for node in dag.op_nodes(op=HGate):
            pre_node = next(dag.quantum_predecessors(node))
            if (pre_node.name == 'h'):
                annihilation()
                fl = False
        # Remove CX Gate
        for node in dag.op_nodes(op = CXGate):
            pre_nodes = list(dag.quantum_predecessors(node))
            if (pre_nodes[0].name == 'cx' and len(pre_nodes) == 1 and pre_nodes[0].qargs == node.qargs):
                dag.remove_op_node(node)
                dag.remove_op_node(pre_nodes[0])
                fl = False
        if fl:
            break
    # dag_drawer(dag)
    return dag_to_circuit(dag)

def read_esop(filepath):
    with open(filepath, "r") as fp:
        lines = fp.readlines()

    n_pi   = 0
    n_po   = 0
    n_cube = 0
    esop_c = []
    esop_x = []

    for line in lines:
        if line[0] == '#':
            continue
        elif line[0] == '.':
            if line[1] == 'i':
                n_pi = int(line[3:])
            if line[1] == 'o':
                n_po = int(line[3:])
            if line[1] == 'p':
                n_cube = int(line[3:])
        else:
            esop_c.append(line.strip().split(' ')[0])
            esop_x.append(line.strip().split(' ')[1])
        
    return {'n_pi':   n_pi,
            'n_po':   n_po,
            'n_cube': n_cube,
            'esop_c': esop_c,
            'esop_x': esop_x}

# Add Expand Toffoli
def add_toffoli(qc, c1, c2, x):
    qc.h(x)
    qc.cx(c2, x)
    qc.tdg(x)
    qc.cx(c1, x)
    qc.t(x)
    qc.cx(c2, x)
    qc.tdg(x)
    qc.cx(c1, x)
    qc.t(c2)
    qc.t(x)
    qc.cx(c1, c2)
    qc.h(x)
    qc.t(c1)
    qc.tdg(c2)
    qc.cx(c1, c2)

# Add Expand Toffoli Reverse
def add_toffoli_r(qc, c1, c2, x):
    qc.cx(c1, c2)
    qc.tdg(c2)
    qc.t(c1)
    qc.h(x)
    qc.cx(c1, c2)
    qc.t(x)
    qc.t(c2)
    qc.cx(c1, x)
    qc.tdg(x)
    qc.cx(c2, x)
    qc.t(x)
    qc.cx(c1, x)
    qc.tdg(x)
    qc.cx(c2, x)
    qc.h(x)

# ESOP -> MCX
def esop_mcx(filepath):
    esop   = read_esop(filepath)
    n_pi   = esop['n_pi']
    n_po   = esop['n_po']
    n_cube = esop['n_cube']
    esop_c = esop['esop_c']
    esop_x = esop['esop_x']

    qr = QuantumRegister(2 * n_pi + n_po - 1)
    qc = QuantumCircuit(qr)

    for i in range(n_cube):
        bits_c = []
        bits_x = []
        esop_x[i]
        for idx in range(n_pi):
            if esop_c[i][idx] == '0':
                qc.x(idx)
            if esop_c[i][idx] == '0' or esop_c[i][idx] == '1':
                bits_c.append(idx)
        for idx in range(n_po):
            if esop_x[i][idx] == '1':
                bits_x.append(idx)
                qc.mcx(bits_c, idx + 2 * n_pi - 1)
        for idx in range(n_pi):
            if esop_c[i][idx] == '0':
                qc.x(idx)   

    return qc

# ESOP -> CCX
def esop_ccx(filepath):
    esop   = read_esop(filepath)
    n_pi   = esop['n_pi']
    n_po   = esop['n_po']
    n_cube = esop['n_cube']
    esop_c = esop['esop_c']
    esop_x = esop['esop_x']

    qr = QuantumRegister(2 * n_pi + n_po - 1)
    qc = QuantumCircuit(qr)

    for i in range(n_cube):
        bits_c  = [i for i, x in enumerate(esop_c[i]) if x == '0' or x == '1']
        bits_cp = [i for i, x in enumerate(esop_c[i]) if x == '1']
        bits_cn = [i for i, x in enumerate(esop_c[i]) if x == '0']
        bits_x  = [i for i, x in enumerate(esop_x[i]) if x == '1']
        
        for idx in bits_cn:
            qc.x(idx)
        
        if len(bits_c) == 1:
            for idx in bits_x:
                qc.cx(bits_c[0], 2 * n_pi + idx - 1)
        elif len(bits_c) == 2:
            for idx in bits_x:
                qc.ccx(bits_c[0], bits_c[1], 2 * n_pi + idx - 1)
        else:
            idx_ancilla = n_pi
            qc.ccx(bits_c[0], bits_c[1], idx_ancilla)
            for idx in bits_c[2:]:
                qc.ccx(idx, idx_ancilla, idx_ancilla + 1)
                idx_ancilla = idx_ancilla + 1
            
            for idx in bits_x:
                qc.cx(idx_ancilla, 2 * n_pi + idx - 1)
                
            for idx in bits_c[:1:-1]:
                qc.ccx(idx, idx_ancilla - 1, idx_ancilla)
                idx_ancilla = idx_ancilla - 1
            qc.ccx(bits_c[0], bits_c[1], idx_ancilla)

        for idx in bits_cn:
            qc.x(idx)

    return qc

# ESOP -> CX, trivial
def esop_cx_trivial(filepath):
    esop   = read_esop(filepath)
    n_pi   = esop['n_pi']
    n_po   = esop['n_po']
    n_cube = esop['n_cube']
    esop_c = esop['esop_c']
    esop_x = esop['esop_x']

    qr = QuantumRegister(2 * n_pi + n_po - 1)
    qc = QuantumCircuit(qr)

    for i in range(n_cube):
        bits_c  = [i for i, x in enumerate(esop_c[i]) if x == '0' or x == '1']
        bits_cp = [i for i, x in enumerate(esop_c[i]) if x == '1']
        bits_cn = [i for i, x in enumerate(esop_c[i]) if x == '0']
        bits_x  = [i for i, x in enumerate(esop_x[i]) if x == '1']
        
        for idx in bits_cn:
            qc.x(idx)
        
        if len(bits_c) == 1:
            for idx in bits_x:
                qc.cx(bits_c[0], 2 * n_pi + idx - 1)
        elif len(bits_c) == 2:
            for idx in bits_x:
                add_toffoli(qc, bits_c[0], bits_c[1], 2 * n_pi + idx - 1)
        elif len(bits_c) > 2:
            idx_ancilla = n_pi
            add_toffoli(qc, bits_c[0], bits_c[1], idx_ancilla)
            for idx in bits_c[2:]:
                add_toffoli(qc, idx, idx_ancilla, idx_ancilla + 1)
                idx_ancilla = idx_ancilla + 1
            
            for idx in bits_x:
                qc.cx(idx_ancilla, 2 * n_pi + idx - 1)
                
            for idx in bits_c[:1:-1]:
                add_toffoli(qc, idx, idx_ancilla - 1, idx_ancilla)
                idx_ancilla = idx_ancilla - 1
            add_toffoli(qc, bits_c[0], bits_c[1], idx_ancilla)

        for idx in bits_cn:
            qc.x(idx)

    return qc

# ESOP -> CX
def esop_cx(filepath):
    esop   = read_esop(filepath)
    n_pi   = esop['n_pi']
    n_po   = esop['n_po']
    n_cube = esop['n_cube']
    esop_c = esop['esop_c']
    esop_x = esop['esop_x']

    qr = QuantumRegister(2 * n_pi + n_po - 1)
    qc = QuantumCircuit(qr)

    for i in range(n_cube):
        bits_c  = [i for i, x in enumerate(esop_c[i]) if x == '0' or x == '1']
        bits_cp = [i for i, x in enumerate(esop_c[i]) if x == '1']
        bits_cn = [i for i, x in enumerate(esop_c[i]) if x == '0']
        bits_x  = [i for i, x in enumerate(esop_x[i]) if x == '1']
        
        for idx in bits_cn:
            qc.x(idx)
        
        if len(bits_c) == 1:
            for idx in bits_x:
                qc.cx(bits_c[0], 2 * n_pi + idx - 1)
        elif len(bits_c) == 2:
            for idx in bits_x:
                add_toffoli(qc, bits_c[0], bits_c[1], 2 * n_pi + idx - 1)
        elif len(bits_c) > 2:
            idx_ancilla = n_pi
            add_toffoli(qc, bits_c[0], bits_c[1], idx_ancilla)
            for idx in bits_c[2:]:
                add_toffoli(qc, idx, idx_ancilla, idx_ancilla + 1)
                idx_ancilla = idx_ancilla + 1
            
            for idx in bits_x:
                qc.cx(idx_ancilla, 2 * n_pi + idx - 1)
                
            for idx in bits_c[:1:-1]:
                add_toffoli_r(qc, idx, idx_ancilla - 1, idx_ancilla)
                idx_ancilla = idx_ancilla - 1
            add_toffoli_r(qc, bits_c[0], bits_c[1], idx_ancilla)

        for idx in bits_cn:
            qc.x(idx)

    return qc

def zx_simp(qc):
    qasm = qc.qasm().replace('pi', '3.14159')
    c      = zx.Circuit.from_qasm(qasm)
    g      = c.to_graph()
    zx.full_reduce(g)
    cc     = zx.extract_circuit(g.copy())
    qasm   = zx.Circuit.to_qasm(cc)
    qc_opt = QuantumCircuit.from_qasm_str(qasm)
    return qc_opt

def op_count(qc):
    return reduce(lambda x, y: x + y[1], list(qc.count_ops().items()), 0)

def test_eq(qc1, qc2):
    uni_sim = Aer.get_backend('unitary_simulator')
    uni1 = execute(qc1, uni_sim).result().get_unitary()
    uni2 = execute(qc2, uni_sim).result().get_unitary()
    # print(uni1 == uni2)
    # print((uni1 - uni2) < 1e-10)
    print(sum(sum((uni1 - uni2) < 1e-10)))
    print(uni1 - uni2)

def run_exp(filepath):
    qc_cxt = esop_cx_trivial(filepath)
    qc_cx  = esop_cx(filepath)
    qc_t   = reduce_t(qc_cx)
    print('HI')
    print(op_count(qc_cxt), op_count(qc_t))
    qc_zxt = zx_simp(qc_cxt)
    print('HI')
    qc_zx  = zx_simp(qc_t)
    print(op_count(qc_cxt), op_count(qc_zxt), op_count(qc_t), op_count(qc_zx))

run_exp('abc/a.esop')
# qc_ccx = esop_ccx('abc/b.esop')
# qc_mcx = esop_mcx('abc/b.esop')
# qc     = esop_cx_trivial('abc/b.esop')
# qc_red = reduce_t(qc)
# test_eq(qc_ref, qc_red)
# test_eq(qc_ref, qc_opt)
# test_eq(qc_ccx, qc_mcx)

# print(qc.draw())
# print(qc_opt.draw())
# print(op_count(qc))
# print(op_count(qc_opt))

