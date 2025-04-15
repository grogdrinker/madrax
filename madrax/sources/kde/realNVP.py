
### from https://github.com/senya-ashukha/real-nvp-pytorch/blob/master/real-nvp-pytorch.ipynb

import torch,math
from torch import nn
import numpy as np

class RealNVP(nn.Module):
	def __init__(self, nfea,device="cpu"):
		super(RealNVP, self).__init__()

		nets = lambda: nn.Sequential(nn.Linear(nfea, 256), nn.LeakyReLU(),
									 nn.Linear(256, 256), nn.LeakyReLU(),
									 nn.Linear(256, nfea), nn.Tanh()).to(device)

		nett = lambda: nn.Sequential(nn.Linear(nfea, 256), nn.LeakyReLU(),
									 nn.Linear(256, 256), nn.LeakyReLU(),
									 nn.Linear(256, nfea)).to(device)
		m = []
		tm = []
		for i in range(nfea):
			if i % 2 == 0:
				tm += [0]
			else:
				tm += [1]
		m += [tm]
		tm = []
		for i in range(nfea):
			if i % 2 == 1:
				tm += [0]
			else:
				tm += [1]
		m += [tm]
		mask = torch.tensor(m * 3).to(device).float()

		self.device = device

		self.mask = nn.Parameter(mask, requires_grad=True)
		self.t = torch.nn.ModuleList([nett() for _ in range(len(mask))])
		self.s = torch.nn.ModuleList([nets() for _ in range(len(mask))])

	def to(self,device):
		self.t = self.t.to(device)
		self.s = self.s.to(device)

		m = self.mask.data.to(device)
		del self.mask
		self.mask = m
		self.device = device

	def g(self, z):
		x = z
		for i in range(len(self.t)):
			x_ = x*self.mask[i]
			s = self.s[i](x_)*(1 - self.mask[i])
			t = self.t[i](x_)*(1 - self.mask[i])
			x = x_ + (1 - self.mask[i]) * (x * torch.exp(s.clamp(max=5)) + t)
		return x

	def f(self, x):
		log_det_J, z = x.new_zeros(x.shape[0],device=self.device), x
		for i in reversed(range(len(self.t))):
			z_ = self.mask[i] * z
			s = self.s[i](z_) * (1-self.mask[i])

			t = self.t[i](z_) * (1-self.mask[i])
			z = (1 - self.mask[i]) * (z - t) * torch.exp(-s.clamp(min=-5)) + z_
			log_det_J -= s.sum(dim=1)
		return z, log_det_J

	def prior_logprob(self,x):
		size = x.shape[1]
		norm_const = 1.0 / (math.pow((2 * math.pi), float(size) / 2))

		modQuad = (x.norm(dim=-1).clamp(max=10)**2).clamp(min=-10)
		return torch.log((torch.exp( -0.5 * modQuad) * norm_const).clamp(min=0.01))

	def log_prob(self,x): #logprob

		z, logp = self.f(x)
		return self.prior_logprob(z) + logp
	def sample(self, batchSize):

		a = torch.randn(batchSize,self.s[0][0].in_features,device=self.device)
		x = self.g(a)
		return x

def train_kde(dataset,epochs = 2001):
	dev = dataset.device
	nfea = dataset.shape[1]

	flow = RealNVP(nfea = nfea,device=dev)
	optimizer = torch.optim.Adam([p for p in flow.parameters() if p.requires_grad==True], lr=1e-4)
	
	for t in range(epochs):
		idx = np.random.randint(0,len(dataset), size=10000)
		
		noisy_moons = dataset[idx]
	
		loss = -flow.log_prob(noisy_moons).mean()

		loss.backward(retain_graph=True)
		optimizer.step()
		optimizer.zero_grad()

		if t % 10 == 0:
			print('iter %s:' % t, 'loss = %.3f' % loss)
			
	return flow
