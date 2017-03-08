int read_sw_lines(const char *data_dir);

void read_space_weather(const char *data_dir, struct msis_struct *pmsis, int nobs, double cf107, double cf107A, double cAp);

int find_msis_time(int nobs, double year, double doy);

void get_msis_ap_array(struct msis_struct *pmsis, int index, double hour, double msis_ap_array[7]);
