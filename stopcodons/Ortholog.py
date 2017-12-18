class Ortholog(object):
	"""docstring for Ortholog"""
	def __init__(self, id1, name1, orthologyType, name2, id2):
		super(Ortholog, self).__init__()
		self.id1 = id1
		self.id2 = id2
		self.name1 = name1
		self.name2 = name2
		self.orthologyType = orthologyType

		