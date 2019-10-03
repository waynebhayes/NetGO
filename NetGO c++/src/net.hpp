
class NetGO{
	public:
		NetGO();
		void newCluster(std::string alignFile);
		void newg2g(std::string g2gFile);
		std::map<std::string, int> SetIntersect(std::map<std::string, int> T1, std::map<std::string, int> T2);
		float K_g(std::string g);
		float K_gset(std::map<std::string, int> T);
		float K_p(std::string p);
		float K_AC(std::map<int, std::vector<std::string>> C);
		void get_pC_CA(std::ifstream& alignFile);
		void get_pGO_GOp(std::ifstream& g2gFile);
		bool DRACONIAN;
		bool VERBOSE;
		std::map<int, std::vector<std::string>> CA;
		std::map<std::string, std::map<int, int>> pC;
		std::map<std::string, std::map<std::string, int>> pGO;
		std::map<std::string, std::map<std::string, int>> GOp;
		std::map<std::string, int> GOfreq;
		std::unordered_set<std::string> encountered_proteins;
};
