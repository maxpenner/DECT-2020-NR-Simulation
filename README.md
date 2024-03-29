# DECT-2020 New Radio Link-Level Simulation Environment
This program contains a Matlab link-level simulation environment for the ETSI standard DECT-2020 New Radio (NR). DECT-2020 NR is a Radio Interface Technology (RIT) designed to provide a slim but powerful technology foundation for wireless applications deployed in various use cases and markets.

[DECT Introduction](https://www.etsi.org/technologies/dect)

The standard consists of multiple parts which can be found on the DECT Technical Committee (TC) webpage.

[DECT committee page](https://www.etsi.org/committee/1394-dect)

## Capabilities
The complete physical layer of a DECT-2020 NR transmitter is implemented. This includes all bandwidths, MIMO modes, channel coding etc. Additionally, BERs and PERs can be simulated in different wireless channel models, in particular doubly-selective channels. For the receiver, STO and CFO synchronization as well as most MIMO modes have been implemented.

## Main Scripts
- **main_single_packet.m**: Simple example demonstrating how to use the simulation environment. Creates a DECT-2020 NR packet, sends it through a wireless channel and decodes it.
- **main_BER_PER_over_MCS.m**: Parallel simulation of packet error rates over SNRs.
- **main_BER_PER_over_MCS_convergence.m**: Same as **main_BER_PER_over_MCS.m**, but optimized for fast convergence at each SNR.
- **main_BER_PER_over_MCS_convergence_full.m**: Same as **main_BER_PER_over_MCS_convergence.m**, but including resampling, an ADC model as well as synchronization of STO and CFO.
- **main_BER_PER_over_MCS_plot_PCC.m**: Plot simulation results of above scripts for PCC.
- **main_BER_PER_over_MCS_plot_PDC.m**: Plot simulation results of above scripts for PDC.

## Requirements
The Matlab LTE Toolbox is required for channel coding, and the Parallel Computing Toolbox for reducing simulation time.

## Exemplary Results

### Packet Error Rates (PERs)
PERs of a SIMO (two receive antennas) system for different MCS over SNR in a Rician fading channel.
<p align="center">
  <img src="./pics/PER_over_SNR.jpg" width="700">
</p>

### Resource Mapping
Resource mapping of STF, DRS, PCC and PDC in the time-frequency lattice.
<p align="center">
  <img src="./pics/Resource_Mapping.jpg" width="700">
</p>

### Channel Interpolation
Interpolated average path gains for doubly selective channel.
<p align="center">
  <img src="./pics/Channel_Interpolation.jpg" width="700">
</p>
