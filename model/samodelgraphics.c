
#include "samodelgraphics.h"

#define DEBUG 0

#if PLOTTING

void samodel_graphics(model_data *md, float pagesize, int linewidth) {

#if PGPLOT

	int k_scene, k_band, k_region;
	float error; 
	double Rrs_max, kd_min, kd_max, perc1, perc2, K_max;

	// Load PGPLOT. 

  	error = cpgopen("?");
  	if (error < 1) {
    	printf("\nERROR: cannot load graphics library.\n");
    	return ;
  	}

	// Set the page size and aspect ratio to 1.0. 

	cpgpap(pagesize, 1.0);

	// Set line width. 

	cpgslw(linewidth);

	// Set text size. 

	cpgsch(2.5);

 	// Background and text colours. 

 	cpgscr(0, 1.0, 1.0, 1.0);
 	cpgscr(1, 0.0, 0.0, 0.0);

 	// Specify number of panels.

 	cpgsubp(2, 4);

 	// 
 	//	Rrs_measured vs Rrs_modelled at pixel 
 	//

	// Set size of plotting window. 

	Rrs_max = 0.0;
	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[0]; k_band++) {
				if (md->Rrs_measured[k_region][k_scene][k_band] > Rrs_max) {
					Rrs_max = md->Rrs_measured[k_region][k_scene][k_band];
				}
				if (md->Rrs_modelled[k_region][k_scene][k_band] > Rrs_max) {
					Rrs_max = md->Rrs_modelled[k_region][k_scene][k_band];
				}
			}
		}
	}

	Rrs_max *= 110.0;

	cpgenv(400.0, 700.0, 0.0, Rrs_max, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "Rrs (100 x sr\\u-1\\d)", "Rrs\\dmeasured\\u v Rrs\\dmodelled\\u");

	cpgbbuf(); // Begin buffering.

	// Rrs_measured

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
			cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->Rrs_measured[md->origin][k_scene][k_band]);
			cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_measured[md->origin][k_scene][k_band + 1]);
			cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_measured[md->origin][k_scene][k_band + 1]);
			// 
			cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->Rrs_measured[md->origin][k_scene][k_band], 16);
		}

		cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->Rrs_measured[md->origin][k_scene][md->n_bands[k_scene] - 1], 16);
	}

	cpgsci(1);

	// Rrs_modelled

	cpgsls(2);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
			cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->Rrs_modelled[md->origin][k_scene][k_band]);
			cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_modelled[md->origin][k_scene][k_band + 1]);
			cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_modelled[md->origin][k_scene][k_band + 1]);
			// 
			cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->Rrs_modelled[md->origin][k_scene][k_band], 16);
		}

		cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->Rrs_modelled[md->origin][k_scene][md->n_bands[k_scene] - 1], 16);
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering. 	


 	// 
 	//	Rrs_measured vs Rrs_modelled at neighbours
 	//

	// Set size of plotting window. 

	cpgenv(400.0, 700.0, 0.0, Rrs_max, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "Rrs (100 x sr\\u-1\\d)", "All Rrs\\dmeasured\\u v Rrs\\dmodelled\\u in Region");

	cpgbbuf(); // Begin buffering.

	// Rrs_measured

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
				cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->Rrs_measured[k_region][k_scene][k_band]);
				cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_measured[k_region][k_scene][k_band + 1]);
				cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_measured[k_region][k_scene][k_band + 1]);
				// 
				cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->Rrs_measured[k_region][k_scene][k_band], 16);
			}

			cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->Rrs_measured[k_region][k_scene][md->n_bands[k_scene] - 1], 16);
		}
	}

	cpgsci(1);

	// Rrs_modelled

	cpgsls(2);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
				cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->Rrs_modelled[k_region][k_scene][k_band]);
				cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_modelled[k_region][k_scene][k_band + 1]);
				cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->Rrs_modelled[k_region][k_scene][k_band + 1]);
				// 
				cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->Rrs_modelled[k_region][k_scene][k_band], 16);
			}

			cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->Rrs_modelled[k_region][k_scene][md->n_bands[k_scene] - 1], 16);
		}
	}

	cpgsls(1);
	cpgsci(1);

	cpgebuf(); // End buffering. 	

	//
	// Percentage relative RMS error measured vs modelled at pixel in each scene.
	//

	cpgenv(400.0, 700.0, 0.0, 100.0, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "RMS Error (%)", "Rrs\\dmeasured\\u v Rrs\\dmodelled\\u - Relative Error (%)");

	cpgbbuf(); // Begin buffering.

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
			perc1 = 100.0*fabs(md->Rrs_modelled[md->origin][k_scene][k_band] - md->Rrs_measured[md->origin][k_scene][k_band])/md->Rrs_measured[md->origin][k_scene][k_band];
			perc2 = 100.0*fabs(md->Rrs_modelled[md->origin][k_scene][k_band + 1] - md->Rrs_measured[md->origin][k_scene][k_band + 1])/md->Rrs_measured[md->origin][k_scene][k_band + 1];
			cpgmove(md->wavelengths[k_scene][k_band],     perc1);
			cpgdraw(md->wavelengths[k_scene][k_band + 1], perc2);
			cpgmove(md->wavelengths[k_scene][k_band + 1], perc2);
			// 
			cpgpt1(md->wavelengths[k_scene][k_band], perc1, 16);
		}

		perc1 = 100.0*fabs(md->Rrs_modelled[md->origin][k_scene][md->n_bands[k_scene] - 1] - md->Rrs_measured[md->origin][k_scene][md->n_bands[k_scene] - 1])/md->Rrs_measured[md->origin][k_scene][md->n_bands[k_scene] - 1];
		cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], perc1, 16);
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering.	


	//
	// Percentage relative RMS error measured vs modelled at neighbours. 
	//

	cpgenv(400.0, 700.0, 0.0, 100.0, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "RMS Error (%)", "All Rrs\\dmeasured\\u v Rrs\\dmodelled\\u in Region - Relative Error (%)");

	cpgbbuf(); // Begin buffering.

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
				perc1 = 100.0*fabs(md->Rrs_modelled[k_region][k_scene][k_band] - md->Rrs_measured[k_region][k_scene][k_band])/md->Rrs_measured[k_region][k_scene][k_band];
				perc2 = 100.0*fabs(md->Rrs_modelled[k_region][k_scene][k_band + 1] - md->Rrs_measured[k_region][k_scene][k_band + 1])/md->Rrs_measured[k_region][k_scene][k_band + 1];
				cpgmove(md->wavelengths[k_scene][k_band],     perc1);
				cpgdraw(md->wavelengths[k_scene][k_band + 1], perc2);
				cpgmove(md->wavelengths[k_scene][k_band + 1], perc2);
				// 
				cpgpt1(md->wavelengths[k_scene][k_band], perc1, 16);
			}

			perc1 = 100.0*fabs(md->Rrs_modelled[k_region][k_scene][md->n_bands[k_scene] - 1] - md->Rrs_measured[k_region][k_scene][md->n_bands[k_scene] - 1])/md->Rrs_measured[k_region][k_scene][md->n_bands[k_scene] - 1];
			cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], perc1, 16);
		}
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering.	

	//
	// Rrs deep 
	//

	// Set size of plotting window. 

	Rrs_max = 0.0;
	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[0]; k_band++) {
				if (md->rrs_dp[k_region][k_scene][k_band] > Rrs_max) {
					Rrs_max = md->rrs_dp[k_region][k_scene][k_band];
				}
			}
		}
	}

	Rrs_max *= 110.0;

	cpgenv(400.0, 700.0, 0.0, Rrs_max, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "rrs (100 x sr\\u-1\\d)", "All rrs\\d\\(0766)\\u in Region");

	cpgbbuf(); // Begin buffering.

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
				cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->rrs_dp[k_region][k_scene][k_band]);
				cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->rrs_dp[k_region][k_scene][k_band + 1]);
				cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->rrs_dp[k_region][k_scene][k_band + 1]);
				// 
				cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->rrs_dp[k_region][k_scene][k_band], 16);
			}

			cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->rrs_dp[k_region][k_scene][md->n_bands[k_scene] - 1], 16);
		}
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering.

	//
	// Rrs bottom 
	//

	// Set size of plotting window. 

	Rrs_max = 0.0;
	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[0]; k_band++) {
				if (md->rrs_bottom[k_region][k_scene][k_band] > Rrs_max) {
					Rrs_max = md->rrs_bottom[k_region][k_scene][k_band];
				}
			}
		}
	}

	Rrs_max *= 110.0;

	cpgenv(400.0, 700.0, 0.0, Rrs_max, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "rrs (100 x sr\\u-1\\d)", "All rrs\\dB\\u in Region");

	cpgbbuf(); // Begin buffering.

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_region = 0; k_region < md->n_regions; k_region++) {
			for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
				cpgmove(md->wavelengths[k_scene][k_band],     100.0*md->rrs_bottom[k_region][k_scene][k_band]);
				cpgdraw(md->wavelengths[k_scene][k_band + 1], 100.0*md->rrs_bottom[k_region][k_scene][k_band + 1]);
				cpgmove(md->wavelengths[k_scene][k_band + 1], 100.0*md->rrs_bottom[k_region][k_scene][k_band + 1]);
				// 
				cpgpt1(md->wavelengths[k_scene][k_band], 100.0*md->rrs_bottom[k_region][k_scene][k_band], 16);
			}

			cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], 100.0*md->rrs_bottom[k_region][k_scene][md->n_bands[k_scene] - 1], 16);
		}
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering.


	//
	// K 
	//

	// Set size of plotting window. 

	K_max = 0.0;
	for (k_region = 0; k_region < md->n_regions; k_region++) {
		for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
			for (k_band = 0; k_band < md->n_bands[0]; k_band++) {
				if (md->K[k_scene][k_band] > K_max) {
					K_max = md->K[k_scene][k_band];
				}
			}
		}
	}

	K_max *= 1.1;

	cpgenv(400.0, 700.0, 0.0, K_max, 0, 0);

	// Plot labels. 

	cpglab("\\(0637) (nm)", "K (m\\u-1\\d)", "All K in Region");

	cpgbbuf(); // Begin buffering.

	cpgsls(1);

	for (k_scene = 0; k_scene < md->n_scenes; k_scene++) {
		cpgsci(2 + 2*k_scene);
		for (k_band = 0; k_band < md->n_bands[k_scene] - 1; k_band++) {
			cpgmove(md->wavelengths[k_scene][k_band],     md->K[k_scene][k_band]);
			cpgdraw(md->wavelengths[k_scene][k_band + 1], md->K[k_scene][k_band + 1]);
			cpgmove(md->wavelengths[k_scene][k_band + 1], md->K[k_scene][k_band + 1]);
			// 
			cpgpt1(md->wavelengths[k_scene][k_band], md->K[k_scene][k_band], 16);
		}

		cpgpt1(md->wavelengths[k_scene][md->n_bands[k_scene] - 1], md->K[k_scene][md->n_bands[k_scene] - 1], 16);
	}

	cpgsci(1);
	cpgsls(1);

	cpgebuf(); // End buffering.


	//
	// Data display
	//

	// Set size of plotting window. 

	cpgsci(0);
	cpgenv(0.0, 10.0, 0.0, 10.0, 0, 0);
	cpgsci(1);

	char display_str[100];

	if (md->converged) {
  		sprintf(display_str, "Model converged (%d iterations)", md->n_iterations);
		cpgtext(.0, 12.0, display_str);
	} else {
		cpgtext(.0, 12.0, "Model diverged");
	}

	cpgtext(.0 + 0.0, 10.0, "E");
	cpgtext(.0 + 0.9, 10.0, "=");
  	sprintf(display_str, "%.2f (%%)", md->model_error);
  	cpgtext(.0 + 1.35, 10.0, display_str);

	cpgtext(.0 + 0.0, 8.0, "H");
	cpgtext(.0 + 0.9, 8.0, "=");
  	sprintf(display_str, "%.2f (m)", md->depth);
  	cpgtext(.0 + 1.35, 8.0, display_str);

	cpgtext(.0 + 0.0, 6.0, "K\\dmin\\u");
	cpgtext(.0 + 0.9, 6.0, "=");
  	sprintf(display_str, "%.3f (m\\u-1\\d)", md->K_min);
  	cpgtext(.0 + 1.35, 6.0, display_str);

	cpgtext(4.25, 10.0, "B\\dsand\\u");
	cpgtext(4.25 + 1.6, 10.0, "=");
  	sprintf(display_str, "%.2f (%%)", md->B_type_percent[0]);
  	cpgtext(4.25 + 2.05, 10.0, display_str);

	cpgtext(4.25, 8.0, "B\\dseagrass\\u");
	cpgtext(4.25 + 1.6, 8.0, "=");
  	sprintf(display_str, "%.2f (%%)", md->B_type_percent[1]);
  	cpgtext(4.25 + 2.05, 8.0, display_str);

	cpgtext(4.25, 6.0, "B\\dcoral\\u");
	cpgtext(4.25 + 1.6, 6.0, "=");
  	sprintf(display_str, "%.2f (%%)", md->B_type_percent[2]);
  	cpgtext(4.25 + 2.05, 6.0, display_str);

  	cpgclos(); // Close the graphics device.
#endif
}

#endif
