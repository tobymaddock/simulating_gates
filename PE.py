import qutip as qt
import numpy as np
import abc as ABC
import matplotlib.pyplot as plt

def H1_coeff(t, args):
    return 1 * np.exp(1j*t)

## PARAMETERS ##
ldp = 0.019 # lamb dicke parameter
trap_freq = 2*np.pi*220000 #trap frequency (Hz)

## PAULI MATRICES ##
X = qt.sigmax()
Y = qt.sigmay()
Z = qt.sigmaz()
I = qt.qeye(2)
spinUp = qt.sigmap()
spinDown = qt.sigmam()

## COLLECTIVE SPIN OPERATORS ##
S_x = qt.tensor(X, I) + qt.tensor(I, X)
S_y = qt.tensor(Y, I) + qt.tensor(I, Y)
S_z = qt.tensor(Z, I) + qt.tensor(I, Z)

## MEET THE FOCK STATES ##
fock_init = qt.basis(20, 0)
a = qt.destroy(20)
a_dag = qt.create(20)
I_fock = qt.qeye(20)

## Position and momentum quadratures ## 
pos_quadrature = 1j * (a - a_dag)
momentum_quadrature = a + a_dag

## BASIS STATES ## 
# for now I will simply model this as a 2 level system
ket_0 = qt.basis(2, 0)
ket_1 = qt.basis(2, 1)

## SOLVER OPTIONS ##
options = qt.solver.Options()
print(options)
options.nsteps = 1e9
options.atol = 1e-15
options.rtol = 1e-20
options.method = 'bdf'

class tqGate:
    """Basic outline for all of the Hamiltonians I'm going to be looking at;
        the collapse operators used will be common to any gate """
    def __init__(self):
        self.hamiltonian
        self.cops = []

    def set_time_intervals(self, times):
        """ Set time intervals to simulate over"""
        self.times = times

    def verify_intended_behaviour():
        pass
    
    def set_dephasing_rates(self, dephasing_rates):
        """ (For lindblad master equations) set the dephasing rate associated with 
            each ion """
        self.dephasing_rates = dephasing_rates
    
    def set_collapse_operators():
        """ Set dephasing and depolarising collapse operators """
        pass
    
    
    def solve_via_mc():
        """ Function allows for solving of the two qubit gate via monte carlo simulations """
        pass
    
    def solve_via_lindblad(self, starting_state, times):
        """ Function allows for solving of the two qubit gate via lindblad operators
        Assumptions:
            - Morkovian noise
            - noise can be modelled in the form of of simple dephasing/dissipation rates"""
        self.generate_hamiltonian(times)
        
        return qt.mesolve(self.hamiltonian, starting_state, times, self.cops, options=options)
    
    def solve_via_br(self, input_state, times, aops,  ):
        """ Function allows for noise to be related to the noise power spectrum of
            environmental noise, and the corresponding master equation is then solved
            Requirements:
                - Aops - a list of operators coupling the environment to the syste,
                using a provided PSD (note that the aops can be made time dependent"""
        
    
class arrazola_gate(tqGate):
    """ Gate as outlined in https://www.nature.com/articles/s42005-023-01243-8, which using 
        spin spin coupling to generate entanglement, and dynamical decoupling to protect
        from dephasing noise
        Requires 2 fields:
            - A continuous X drive
            - A continuous Y drive
    """
    def __init__(self):
        self.generate_hamiltonian()
        self.compute_tgate()
        self.cops = []

    def generate_hamiltonian(self, times = [0,1,2]):
        """ 
        Generate a hamiltonian given the experimental parameters
        # note # this hamiltonian uses the RWA
        """
        #[H_1, lambda t, args: np.exp(-1j*t)]
        f_k = 5 # need to compute this
        s_s_coupling = 0.5

        self.det = 2 * ldp * trap_freq * ((f_k**2 + 4*ldp**2 * s_s_coupling**2)**(1/2) + 2 * ldp * s_s_coupling)
        
        # states are of the form, s1@s1@sharedmode
        H_1 = 1/2 * ldp * trap_freq * f_k * (qt.tensor(a, S_z))
        H_2 = - 1/2 * ldp**2 * trap_freq * s_s_coupling * qt.tensor((S_z * S_z), I_fock)
        H = [[H_1, lambda t, args: 1 * np.exp(-self.det*1j*t)], H_2]
        self.hamiltonian = H
        
    def compute_tgate(self):
        """
        Compute the gate time (time required for phase pi/8)

        """
        self.t_gate = 2 * np.pi / self.det
        
    def compute_gate_time(self):
        ""
        pass
        
        
    def compute_spin_spin_coupling():
        pass
             
    def set_x_drive_intensity():
        pass
    
    def set_y_drive_intensity():
        pass
    
    

  
class MS_gate(tqGate):
    """ Generic Molmer Sorensen Gate """
    def __init__(self, times):
        self.times = times
        self.cops = [0 * qt.tensor(I, qt.sigmaz(), I_fock), # dephasing dissipation operators
                     0 * qt.tensor(qt.sigmaz(), I, I_fock),
                     0 * qt.tensor(qt.sigmax(), I,  I_fock), # depolarisation dissipation
                     0 *qt.tensor(I, qt.sigmax(), I_fock),
                     0 * qt.tensor(I, I, a),
                     0 * qt.tensor(I, I, a_dag)] # heating
    
        self.ldp = 0.019
        
    def set_dephasing_rates(self, dephasing_rates):
        """ adjust collapse operators according to entered dephasing rates """
        # for now will just assume that the dephasing rate experienced by each ion is the same
        coeff = np.sqrt(dephasing_rates/2)
        self.cops = [coeff * qt.tensor(I, qt.sigmaz(), I_fock), # dephasing dissipation operators
                     coeff * qt.tensor(qt.sigmaz(), I, I_fock)]
    
    def set_heating_rates(self, heating_rates):
        """ adjust collapse operators according to entered dephasing rates """
        # for now will just assume that the dephasing rate experienced by each ion is the same
        coeff = np.sqrt(heating_rates)
        self.cops = [[coeff * qt.tensor(I, I, a),
                     coeff * qt.tensor(I, I, a_dag)]]
    
    
    ## SPIN OPERATOR ##
    def spin_op(self, spin_phase):
        """ Define spin operator for Molmer Sorensen Gate operations """
        return 1j* (spinUp - spinDown)
    
    def generate_hamiltonian(self, times):
        """ Generate a hamiltonian given the experimental parameters
            Assume that ions are in the COM mode, so equal ldps"""
        hbar = 1 # use natural units
        rabi = 60e3
        det  = 2 * self.ldp * rabi
        print(det)
        H_1 = 1j *hbar* self.ldp * rabi/2 * self.spin_op(0)
        H_1 = qt.tensor(H_1, I)
        
        H_2 = 1j *hbar* self.ldp * rabi/2 * self.spin_op(0)
        H_2 = qt.tensor(I, H_2)
        
        H = [[-qt.tensor(H_1, a_dag), lambda t, args: np.exp(1j *det *t)],
             [qt.tensor(H_1, a), lambda t, args: np.exp(-1j* det *t)], 
             [-qt.tensor(H_2,  a_dag), lambda t, args:  np.exp((1j* det *t))],
             [qt.tensor(H_2,  a), lambda t, args: np.exp(-1j* det *t)]]
        
        self.hamiltonian = H

        
        
    def set_experimental_params(self):
        pass
    
    def set_ldp(self, ldp):
        """ Set the assumed Lamb-Dicke parameter"""
        self.ldp = ldp
    
    def compute_cops(self):
        """ Idea is that this function will convert measured (or approximated) noise
            PSDs into a bunch of dephasing/damping collapse operators """
        pass



def infidelity_vs_dephasing():
    """ Plot gate fidelity vs dephasing rates"""
    rabi_freq = 60e3 # chosen rabi frequency
    detuning = 2 * ldp * rabi_freq
    time = (2*np.pi/detuning)
    final_state = 1/np.sqrt(2) *\
        (qt.tensor(ket_0, ket_1) + 1j *qt.tensor(ket_1, ket_0))
    final_state = qt.tensor(final_state, fock_init)

    # phase_acquired = np.pi * ldp**2  * rabi_freq**2 / (2* detuning**2)
    simulation_times = np.linspace(0, time, 2)
    input_state = qt.tensor(ket_0, ket_1, fock_init)

    input_state_as_dm = qt.ket2dm(input_state)
    
    dephasing_rates = np.linspace(0, 3, 100)
    fidelities = []
    for entry in dephasing_rates:
        ms_gate = MS_gate(simulation_times)
        ms_gate.set_dephasing_rates(entry)
        a = ms_gate.solve_via_lindblad(input_state, simulation_times)
        states = a.states
        norm_state = a.states[-1]
        overlap = qt.Qobj.overlap(final_state, norm_state)
        fidelity = np.conjugate(overlap) * overlap
        fidelities.append(fidelity)
        print(fidelity * 100,' ', entry)
        
    return dephasing_rates, fidelities

def infidelity_vs_heating():
    """ Plot gate fidelity vs dephasing rates"""
    rabi_freq = 60e3 # chosen rabi frequency
    detuning = 2 * ldp * rabi_freq
    time = (2*np.pi/detuning)
    final_state = 1/np.sqrt(2) *\
        (qt.tensor(ket_0, ket_1) + 1j *qt.tensor(ket_1, ket_0))
    final_state = qt.tensor(final_state, fock_init)

    # phase_acquired = np.pi * ldp**2  * rabi_freq**2 / (2* detuning**2)
    simulation_times = np.linspace(0, time, 2)
    input_state = qt.tensor(ket_0, ket_1, fock_init)

    input_state_as_dm = qt.ket2dm(input_state)
    
    heating_rates = np.linspace(0, 10, 20)
    fidelities = []
    for entry in heating_rates:
        ms_gate = MS_gate(simulation_times)
        ms_gate.set_heating_rates(entry)
        a = ms_gate.solve_via_lindblad(input_state, simulation_times)
        states = a.states
        norm_state = a.states[-1]
        overlap = qt.Qobj.overlap(final_state, norm_state)
        fidelity = np.conjugate(overlap) * overlap
        fidelities.append(fidelity)
        print(fidelity * 100,' ', entry)
        
    return heating_rates, fidelities


x, y = infidelity_vs_heating()
plt.plot(x, y)

plt.xlabel('COM mode heating rate (phonons/s')
plt.ylabel('Gate Fidelity')
plt.title('Gate Fidelity vs heating rate (COM mode)')

rabi_freq = 60e3 # chosen rabi frequency
detuning = 2 * ldp * rabi_freq
time = (2*np.pi/detuning)
final_state = 1/np.sqrt(2) *\
    (qt.tensor(ket_0, ket_1) + 1j *qt.tensor(ket_1, ket_0))
final_state = qt.tensor(final_state, fock_init)

 # phase_acquired = np.pi * ldp**2  * rabi_freq**2 / (2* detuning**2)
simulation_times = np.linspace(0, time, 5000)
input_state = qt.tensor(ket_0, ket_1, fock_init)

rabi_freq = 60e3 # chosen rabi frequency
ms_gate = MS_gate(simulation_times)
a = ms_gate.solve_via_lindblad(input_state, simulation_times)
states = a.states
norm_state = a.states[-1]
overlap = qt.Qobj.overlap(final_state, norm_state)
fidelity = np.conjugate(overlap) * overlap
print("fidelity is:", fidelity)

# def calc_ldp_for_ion_height(height):
    
#     k = 1.4333e-3
#     grad = 6.77301 * height**(-1.78) * 1.6827*10**(-6)
    
#     return k*grad

# gate_times = []
# ion_heights = []
# fidelities = []

# for i in range(1, 20):   
#     print('ldp' , calc_ldp_for_ion_height(10*i * 1e-6))
#     detuning2 = 2 * calc_ldp_for_ion_height(10*i * 1e-6) * rabi_freq
#     gate_time = (4*np.pi/detuning2)
#     print('gate time (s)', gate_time)
#     gate_times.append(gate_time)
#     ion_heights.append(10*i)
  
    
# plt.plot(ion_heights, gate_times)
# plt.xlabel('Ion Height (um)')
# plt.ylabel('Gate time (s)')
# plt.title('MS Gate Duration vs Ion Height (fixed current)')
# plt.show()


    
final_state = 1/np.sqrt(2) *\
    (qt.tensor(ket_0, ket_1) + 1j *qt.tensor(ket_1, ket_0))
final_state = qt.tensor(final_state, fock_init)


# for entry in ion_heights:
#     print('ldp' , calc_ldp_for_ion_height(10*entry * 1e-6))
#     detuning2 = 2 * calc_ldp_for_ion_height(10*entry * 1e-6) * rabi_freq
#     gate_time = (2*np.pi/detuning2)
        
#     simulation_times = np.linspace(0, gate_time, 200)
#     input_state = qt.tensor(ket_0, ket_1, fock_init)
        
#     input_state_as_dm = qt.ket2dm(input_state)
#     ms_gate = MS_gate(simulation_times)
#     ms_gate.set_ldp(calc_ldp_for_ion_height(10*entry * 1e-6))
#     a_states = ms_gate.solve_via_lindblad(input_state, simulation_times)
        
#     states = a_states.states
#     norm_state = a_states.states[-1]
#     overlap = qt.Qobj.overlap(final_state, norm_state)
#     fidelity = np.conjugate(overlap) * overlap
#     print('fidelity is ', fidelity)
#     fidelities.append(fidelity)
    
# plt.plot(ion_heights, fidelities)
# plt.xlabel('Ion Height (um)')
# plt.ylabel('Gate Fidelity')
# plt.title('MS Gate Fidelity vs Ion Height (fixed current)')
# plt.show()



def generate_molmer_sorensen_graph(states, times):
    upup = qt.tensor(ket_1, ket_1)
    downdown = qt.tensor(ket_0, ket_0)
    updown = qt.tensor(ket_1, ket_0)
    downup = qt.tensor(ket_0, ket_1)
    
    upup_probs = []
    downdown_probs = []
    updown_probs = []
    downup_probs = []
    
    def measurement_prob(input_state, measurement_state = updown):
        proj_upp = qt.tensor(qt.ket2dm(measurement_state), I_fock)* input_state
        measurement_probability = qt.Qobj.overlap(proj_upp, proj_upp)
        return measurement_probability
    
    for i in range(len(states)):
        # normalise state first, then calculate measurement probability
        norm_state = a.states[i]
        #proj_upp = qt.tensor(qt.ket2dm(updown), I_fock)* norm_state
        #measurement_probability = qt.Qobj.overlap(proj_upp, proj_upp)
        updown_probs.append(measurement_prob(norm_state, updown))    
        downup_probs.append(measurement_prob(norm_state, downup))
        upup_probs.append(measurement_prob(norm_state, upup))    
        downdown_probs.append(measurement_prob(norm_state, downdown))
        
    fig, ax = plt.subplots()

    # Plot the first set of data
    ax.plot(simulation_times, updown_probs, label='updown')
    ax.plot(simulation_times, downup_probs, label='downup')
    ax.plot(simulation_times, upup_probs, label='upup')
    ax.plot(simulation_times, downdown_probs, label='downdown')
    ax.set_title('State probabilities over time')
    ax.set_xlabel('time')
    ax.set_ylabel('probability')
    ax.legend()

    # Show the plot
    plt.show()
    return 0

def generate_PST(states, times):
    """ Generate a phase space trajectory over the course of the gate"""
    p = qt.tensor(I, I, pos_quadrature)
    q = qt.tensor(I, I, momentum_quadrature)
    p_expect = []
    q_expect = []
    
    
    for i in range(len(states)):
        # normalise state first, then calculate measurement probability
        norm_state = a.states[i]
        p_expect.append(qt.Qobj.overlap(norm_state, p*norm_state))
        q_expect.append(qt.Qobj.overlap(norm_state, q*norm_state))
        
    fig, ax = plt.subplots()

    # Plot the first set of data
    ax.plot(simulation_times, p_expect, label='p_probs')
   # ax.plot(simulation_times, q_probs, label='q_probs')
    ax.set_title('State probabilities over time')
    ax.set_xlabel('time')
    ax.set_ylabel('probability')
    ax.legend()

    # Show the plot
    plt.show()
    
    pass

generate_molmer_sorensen_graph(states, simulation_times)
generate_PST(states, simulation_times)

#times = np.linspace(0, time, 2000)
def plotPST_noiseless(times):
    ## DISPLACEMENT ## 
    a = [(1j * ldp * rabi_freq/ (2 * detuning)) * np.exp(np.pi/2*1j)*(np.exp(1j  *  detuning * entry) - 1) for entry in times]
    #print(a)
    
    print(rabi_freq* ldp/detuning)
    real_a = [2*np.real(entry) for entry in a]
    print(real_a)
    imag_a = [2*np.imag(entry) for entry in a]
    fig, ax = plt.subplots()
    # Plot the first set of data
    ax.plot(imag_a, real_a)
    ax.set_title('PST without noise')
    ax.set_xlabel('Re(a)')
    ax.set_ylabel('Im(a)')

    # Show the plot
    plt.show()
    return 0 

# plotPST_noiseless(times)

def plotPST_modulated(times):
    ## DISPLACEMENT ## 
    a = [(1j * ldp * rabi_freq/ (2 * detuning)) * np.exp(np.pi/2*1j)*(np.exp(1j  *  detuning * entry) - 1) for entry in times]
    #print(a)
    
    print(rabi_freq* ldp/detuning)
    real_a = [2*np.real(entry) for entry in a]
    print(real_a)
    imag_a = [2*np.imag(entry) for entry in a]
    fig, ax = plt.subplots()
    # Plot the first set of data
    ax.plot(imag_a, real_a)
    ax.set_title('PST without noise')
    ax.set_xlabel('Re(a)')
    ax.set_ylabel('Im(a)')

    # Show the plot
    plt.show()
    return 0 


## Input state is of form |0>|0>|0> ##
#plus_state = (1/np.sqrt(2))*(ket_0 + ket_1)
#minus_state = (1/np.sqrt(2))*(qt.basis(2,0) - qt.basis(2,1))
#input_state = qt.tensor(cat_state, cat_state, fock_init)
#final_state = (1/np.sqrt(2))*(qt.tensor(plus_state, plus_state) 
                             #+ 1j*qt.tensor(minus_state, minus_state))
#final_state = qt.tensor(final_state, fock_init)
#New_Gate = JJ_coupling_gate(simulation_times)

#a = New_Gate.solve_via_lindblad(input_state, simulation_times)
#print(input_state)
#print(final_state)

#for state in a.states:
    #overlap = qt.Qobj.overlap(final_state, a.states[-1])
    #fidelity = np.conjugate(overlap) * overlap
    #print(fidelity)
