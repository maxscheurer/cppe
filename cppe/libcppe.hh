namespace libcppe {

  class CPPE {
    public:
      CPPE();
      ~CPPE();
      void call_pe_energy(const double *densmat);
      bool gen1int_initialized() { return m_gen1int_initialized; }

    private:
      bool m_gen1int_initialized;
      // TODO: build a struct/class with the basis set info?
      int m_nbas;  // number of basis functions
      int m_nnbas; // nbas*(nbas+1)/2
  };

}
