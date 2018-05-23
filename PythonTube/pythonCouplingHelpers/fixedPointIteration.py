import abc
import numpy as np

class IterationScheme:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def iterate(self, x, x_tilde):
        pass


class UnderrelaxationScheme(IterationScheme):

    def __init__(self, underrelaxation_factor):
        self.omega = underrelaxation_factor

    def iterate(self, x, x_tilde):
        res = x_tilde - x
        x = x + self.omega * res
        e = np.linalg.norm(res) / res.size
        return x, e


class IQNILSScheme(IterationScheme):

    col_keep = 20

    def __init__(self, underrelaxation_factor):
        self.omega = underrelaxation_factor
        self.k = 0
        self.res = []
        self.x_tilde = []

    def clear(self, underrelaxation_factor):
        if False: #self.k > 0: # todo usually one should be able to reuse results from last timestep for IQN
            self.omega = underrelaxation_factor
            self.k = 1
            self.res = [self.res[-1]]
            self.x_tilde = [self.x_tilde[-1]]
        else:
            self.omega = underrelaxation_factor
            self.k = 0
            self.res = []
            self.x_tilde = []

    def iterate(self, x, x_tilde):
        k = self.k
        dim = x.size
        if k is 0:  # do underrelaxation in very first iteration
            self.x_tilde.append(x_tilde)
            self.res.append(self.x_tilde[k] - x)
            x = x + self.omega * self.res[k]
            e = np.linalg.norm(self.res[k]) / self.res[k].size
            self.k += 1
            return x, e
        else:  # do IQN-ILS else
            self.x_tilde.append(x_tilde)
            self.res.append(self.x_tilde[k] - x)
            e = np.linalg.norm(self.res[k]) / self.res[k].size
            V = np.zeros([self.res[k].shape[0],min([self.col_keep,dim,k])])
            W = np.zeros([self.res[k].shape[0],min([self.col_keep,dim,k])])
            jj = 0

            for kk in range(max([k-self.col_keep, k-dim,0]),k):  # maximum number of columns in V, W is dimension of the problem!
                V[:,jj] = self.res[kk+1] - self.res[kk]
                W[:,jj] = self.x_tilde[kk+1] - self.x_tilde[kk]
                jj += 1

            Q, U = np.linalg.qr(V)
            lhs = -(Q.T).dot(self.res[k])
            alpha = np.linalg.solve(U, lhs)  # todo what to do if a linalg error appears here?
            dx = W.dot(alpha)
            x = self.x_tilde[k] + dx
            self.k += 1
            return x,e


