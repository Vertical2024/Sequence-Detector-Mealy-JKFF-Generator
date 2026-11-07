from math import ceil, log2
from sympy import symbols
from sympy.logic.boolalg import SOPform

# ==========================================
# 1. 用户设置
# ==========================================
PATTERN = "0100"          
NON_OVERLAPPING = True    # True: 不重叠检测；False: 重叠检测

# ==========================================
# 2. 构造 Moore 状态机
# ==========================================
def next_state(pattern, current_match, bit):
    if NON_OVERLAPPING and current_match == len(pattern):
        current_match = 0  # 不重叠模式：触发后匹配进度清零

    s = pattern[:current_match] + bit
    for k in range(len(pattern), -1, -1):
        if s.endswith(pattern[:k]):
            return k
    return 0

states = list(range(len(PATTERN) + 1))
transitions = {}

for state in states:
    for bit in ['0', '1']:
        ns = next_state(PATTERN, state, bit)
        transitions[(state, bit)] = ns

# ==========================================
# 3. 状态编码与表格化提取
# ==========================================
num_states = len(states)
nff = ceil(log2(num_states))

print(f"检测序列: {PATTERN} (不重叠模式: {NON_OVERLAPPING})", flush=True)
print(f"需要触发器数量 (Flip-Flops Needed) = {nff}", flush=True)
print("\n状态编码 (State Encoding)", flush=True)

state_codes = {}
for s in states:
    state_codes[s] = format(s, f'0{nff}b')
    print(f"S{s} = {state_codes[s]}", flush=True)

# JK 激励表
JK = {
    ('0', '0'): ('0', 'X'),
    ('0', '1'): ('1', 'X'),
    ('1', '0'): ('X', '1'),
    ('1', '1'): ('X', '0')
}

vars_list = [symbols(f'Q{i}') for i in range(nff)]
X = symbols('X')
all_vars = vars_list + [X]

J_minterms = [[] for _ in range(nff)]
J_dc       = [[] for _ in range(nff)]
K_minterms = [[] for _ in range(nff)]
K_dc       = [[] for _ in range(nff)]

for state in states:
    current = state_codes[state]
    for bit in [0, 1]:
        ns = transitions[(state, str(bit))]
        nextcode = state_codes[ns]
        row = [int(x) for x in current] + [bit]
        idx = int(''.join(map(str, row)), 2)

        for i in range(nff):
            q = current[i]
            qn = nextcode[i]
            J, K = JK[(q, qn)]

            if J == '1': J_minterms[i].append(idx)
            elif J == 'X': J_dc[i].append(idx)

            if K == '1': K_minterms[i].append(idx)
            elif K == 'X': K_dc[i].append(idx)

# 加入未使用的状态作为 Don't Care (无关项)
total_inputs = 2**(nff + 1)
used_states = set()
for s in states:
    code = state_codes[s]
    for bit in [0, 1]:
        used_states.add(int(code + str(bit), 2))

unused = set(range(total_inputs)) - used_states

for i in range(nff):
    J_dc[i] = list(set(J_dc[i] + list(unused)))
    K_dc[i] = list(set(K_dc[i] + list(unused)))
    J_minterms[i] = list(set(J_minterms[i]))
    K_minterms[i] = list(set(K_minterms[i]))

# ==========================================
# 4. 化简并输出方程
# ==========================================
print("\n==============================", flush=True)
print("JK EQUATIONS", flush=True)
print("==============================", flush=True)

Jexprs = []
Kexprs = []

for i in range(nff):
    Jexpr = SOPform(all_vars, J_minterms[i], J_dc[i])
    Kexpr = SOPform(all_vars, K_minterms[i], K_dc[i])
    Jexprs.append(Jexpr)
    Kexprs.append(Kexpr)
    print(f"\nJ{i} = {Jexpr}", flush=True)
    print(f"K{i} = {Kexpr}", flush=True)

# 输出方程 Z
final_state = state_codes[len(PATTERN)]
z_terms = [f"Q{i}" if b == '1' else f"~Q{i}" for i, b in enumerate(final_state)]

print("\n==============================", flush=True)
print("OUTPUT", flush=True)
print("==============================", flush=True)
print("Z =", " & ".join(z_terms), flush=True)

# ==========================================
# 5. 动态自动化验证 (Dynamic Auto-Verification)
# ==========================================
print("\n==============================", flush=True)
print("VERIFICATION RESULT", flush=True)
print("==============================", flush=True)

all_passed = True

for s in states:
    current_code = state_codes[s]
    q_vals = [int(b) for b in current_code]

    for bit in [0, 1]:
        expected_ns = transitions[(s, str(bit))]
        expected_code = state_codes[expected_ns]

        # 构造 SymPy 变量替换字典
        sub_dict = {vars_list[i]: bool(q_vals[i]) for i in range(nff)}
        sub_dict[X] = bool(bit)

        actual_q_next = []
        for i in range(nff):
            # 动态代入生成的 J, K 方程
            j_val = bool(Jexprs[i].subs(sub_dict))
            k_val = bool(Kexprs[i].subs(sub_dict))
            q_curr = q_vals[i]

            # 根据 JK 特征方程计算次态: Q+ = J*~Q + ~K*Q
            q_next = int((j_val and not q_curr) or (not k_val and q_curr))
            actual_q_next.append(str(q_next))

        actual_code = "".join(actual_q_next)

        if actual_code != expected_code:
            print(f"[错误] 状态 S{s}({current_code}) 在输入 X={bit} 时不符合预期!", flush=True)
            print(f"       期望次态: {expected_code}, 方程计算次态: {actual_code}", flush=True)
            all_passed = False

if all_passed:
    print("验证通过", flush=True)