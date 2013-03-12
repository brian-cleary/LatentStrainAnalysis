class LSA(object):
	def __init__(self,input_path,output_path):
		super(LSA,self).__init__()
		self.input_path = input_path
		self.output_path = output_path
		self.hpfx = 'k, bins: ['