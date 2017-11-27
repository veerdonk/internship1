class dNdSEntry(object):
	'''Object describing the dnds entries'''
	def __init__(self, ref, ort, dnds, dn, ds, name="None"):
		super(dNdSEntry, self).__init__()
		self.ref = ref
		self.ort = ort
		self.dnds = dnds
		self.dn = dn
		self.ds = ds
		self.name = name

	def __str__(self):
		return ("{}\t{}\t{}\t{}\t{}\n"\
			.format(self.ref, self.ort, self.dnds, self.dn, self.ds))

	def __eq__(self, other):
		'''Equals'''
		if other is None:
			return False
		return self.ref == other.ref

	def __key(self):
		'''unique key function used by hash'''
		return (self.ref) #, self.ort, self.dnds

	def __hash__(self):
		'''Hash function'''
		return hash(self.__key())