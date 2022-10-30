import numpy as np

def cubically_equivalent_quaternions (q1,q2):
    
    def q_mult(q1, q2):
        """
        
        
        Parameters
        ----------
        q1, q1: numpy ndarray(4)
            Quaternions. First item is the scalar part
        
        Returns
        -------
        q : numpy ndarray(4)
            Quaternion product
        """
        q = np.zeros(4)
        q[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3]
        q[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2]
        q[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1]
        q[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]
        return (q)

    
    
    q=q_mult(q1, q2)
    
    a0,a1,a2,a3=q[0],q[1],q[2],q[3]
    
    cubic_sym=np.array([
    np.array([a0, a1, a2, a3]),
    np.array([a1, a0, a3,a2]),
    np.array([a2,- a3, a0, a1]),
    np.array([a3, a2, - a1, a0]),
    0.5*np.array([a0 - a1 - a2 - a3, a0 + a1 + a2 - a3, a0 - a1 + a2 + a3, a0 + a1 - a2 + a3]) ,
    0.5*np.array([a0 + a1 + a2 + a3,-a0 + a1 - a2 + a3, - a0 + a1 + a2 - a3, - a0 - a1 + a2 + a3]),
    0.5*np.array([a0 - a1 + a2 - a3, a0 + a1 + a2 + a3, -a0 - a1 + a2 + a3, a0 - a1 - a2 + a3]) ,
    0.5*np.array([a0 + a1 - a2 + a3,-a0 + a1 - a2 - a3, a0 + a1 + a2 - a3, - a0 + a1 + a2 + a3]) ,
    0.5*np.array([a0 + a1 - a2 - a3,-a0 + a1 + a2 - a3, a0 - a1 + a2 - a3, a0 + a1 + a2 + a3] ),
    0.5*np.array([a0 - a1 + a2 + a3, a0 + a1 - a2 + a3,-a0 + a1 + a2 + a3, - a0 - a1 - a2 + a3]),
    0.5*np.array([a0 + a1 + a2 - a3,-a0 + a1 + a2 + a3,-a0 - a1 + a2 - a3, a0 - a1 + a2 + a3]),
    0.5*np.array([a0 - a1 - a2 + a3, a0 + a1 - a2 - a3, a0 + a1 + a2 + a3, - a0 + a1 - a2 + a3]),
    (1/np.sqrt(2))*np.array([a0 - a1, a0 + a1, a2 + a3,-a2 + a3]) ,
    (1/np.sqrt(2))*np.array([a0 - a2, a1 - a3, a0 + a2, a1 + a3]) ,
    (1/np.sqrt(2))*np.array([a0 - a3, a1 + a2, -a1 + a2, a0 + a3]) ,
    (1/np.sqrt(2))*np.array([-a1 - a2, a0 - a3, a0 + a3, a1 - a2]),
    (1/np.sqrt(2))*np.array([-a2 - a3, a2 - a3, a0 - a1, a0 + a1]) ,
    (1/np.sqrt(2))*np.array([-a1 - a3, a0 + a2, -a1 + a3, a0 - a2]),
    (1/np.sqrt(2))*np.array([a0 + a1, -a0 + a1, a2 - a3, a2 + a3]),
    (1/np.sqrt(2))*np.array([a0 + a2, a1 + a3, -a0 + a2,- a1 + a3]),
    (1/np.sqrt(2))*np.array([a0 + a3, a1 - a2, a1 + a2, - a0 + a3]) ,
    (1/np.sqrt(2))*np.array([a1 - a2, - a0 - a3, a0 - a3, a1 + a2] ),
    (1/np.sqrt(2))* np.array([a2 - a3, a2 + a3, - a0 - a1, a0 - a1] ) ,
    (1/np.sqrt(2))*np.array([a1 - a3, - a0 + a2, - a1 - a3, a0 + a2] )
    
    ])
    cos_half_theta_max = 0.0
    for i in range(len(cubic_sym)):
        cos_half_theta = np.abs(q_mult(cubic_sym[i], q1).dot(q2))
        if cos_half_theta > cos_half_theta_max:
                    cos_half_theta_max = cos_half_theta
    if cos_half_theta_max > 1.0:
                mis = 0.0
    else:
        mis = 2 * np.arccos(cos_half_theta_max)
        
        
        
    return(mis)