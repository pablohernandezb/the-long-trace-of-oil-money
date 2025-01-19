## The Long Trace of Oil Money Project

This repository contains the code and data for my research project on **a social network analysis on the Chavismo.inc dataset created by Transparencia Venezuela and others**.

**1. Introduction**

* **Research Question:** Who are the central actors in this corruption network? How resilient is this network? And on a more normative dimension, what are the policy recommendations for reforms to be implemented in institutions in order to prevent for incurring again in corruption?
* **Objective:** The object of this study is two-fold with both academic and political implications. First, on the academic side, I expect to test theoretical expectations of social network analysis (SNA) using the case of corruption networks in Venezuela, where I expect to identify actors that are important in the networks based on: (i) The importance on the structure if they are leaders and help organize and direct the corruption; (ii) If they are part of communication channels based on their location on the network; and (iii) assess the likelihood of being prosecuted internationally by the type of the types formed and a given actor's characteristics.
* **Significance:** Following the approach by Cunningham et al. (2016) in Understanding Dark Networks and applying to the Venezuelan case this study aims to shed lights on the way this networks operates and more importantly evaluate the way in which actors in the international community has been dealing with this issue.
* **Background:** Previous studies are not clear in determining the actual brokers of the network and delving more into its functioning. Another important contribution is to assess the international community's behavior when engages in disruption of corruption network vis à vis the US' behavior and point out any similarities and differences.

**2. Data**

* **Data Sources:**
    * The data expected to use in this study to start solving this puzzle will come mainly from Transparency International (2020) and the report called Chavismo Inc. made in an alliance with other NGO's, Alianza Rebelde Investiga (ARI) and the journalism Latin American platform CONNECTAS. 
* **Data Description:**
    * This report includes a social network that covers 86 cases of investigation in 61 countries and encompass 751 agents (persons of interest), 239 institutions. Taken together, this social network includes over 3,900 relationships that can vary and range from:
        * Occupied functions either in the private or public sector
        * Corruption denunciation, by either an NGO or a governmental agency
        * Enablers, helped the main agents in their activity
        * Family, showing kinship between actors
        * Business connection, whether it is part of a company or have a business with another company
        * International trials, show the active cases of actors
        * Company connections, if the agent has dealings with other companies
        * Sanctioned by, agents that received any type of sanction by a country
        * Contract, show agents that subscribed contracts for services
        * Friends, connections on the social side
        * Student, if the actors know each other from college or high school
        * Enemies, shows if the actors hold a grudge
        * Integrates company, a tie when the actor works within a given company
        * Human rights violation, shows a tie when an actor received a claim for violations of human rights;
        * Designates, when a public officer designates an actor to perform public service in either a government agency or a state company 

**3. Methodology**

* **Research Design:** 
    * First, an overall characterization of the network and the construction of basic and topographic network-level measures such as centralization. Moreover, another actor-level measurements will be included to provide further knowledge on the prestige and role a given member fulfills in the network. 
    * A comprehensive explanation on the determinants for tie formation in the network based on the Quadratic Assignment Procedure (QAP) results and logistic regression on the possibility of having either an International or US trial for corruption.
    * A logistics regressions are performed at the actor-level to determine the differences and similarities on the strategies taken by both International and US authorities in order to diminish the operations and effects of the corruption network that originated in Venezuela.
* **Software and Packages:**
    * In this project I employed R, RStudio, tidyverse, ggplot2, statnet, Stata, Excel mainly for data cleaning and merging.
    * For the network analysis I will use the Statnet package developed in R and the associated libraries that have a wide range of functions that allows to characterize and describe a network.

**4. Code**

* **Structure:** 
    * The repository contains the following scripts: 
        * `main.R`: Performs the main social network analyses.
        * `main.do`: Creates variables and some exploratory logistic regression models to understand corruption dynamics.
        
**5. Results**

* **Key Findings:** 
    * [Summarize the main findings of the analysis, e.g., "The results of the logistic regression analysis showed a statistically significant positive relationship between social media use and voter turnout among young adults."]
    * [Discuss the implications of the results, e.g., "These findings suggest that social media can play a valuable role in mobilizing young voters and increasing political participation."]
* **Visualizations:** 
    * [If applicable, include links to visualizations or plots generated by the code, e.g., "See the following plots for visualizations of the key findings: 
        * [Link to plot 1]
        * [Link to plot 2]"]

**6. Conclusion**

* This study represents a departure from previous studies of corruption in Venezuela in two aspects:
    * First, the model includes more measurements for actor's prestige within the network like cutpoint strength and brokerage level. Second, the behavior of law enforcement agency from the international community and the United States when engage in dark networks disruption and which strategy they used. The way in which prestige is measure using cutpoint strength and brokerage level on top of centrality, offered a more accurate approach to signal who are the most important actors within the network.
    * Moreover, cutpoint strength and brokerage level signal who are the important actors in terms of organization and operation, and are more likely to be targeted for law enforcement internationally and by the US in order to diminish the corrupt activity of such a network.

**7. Contributing**

* Contributions to this project are welcome. You can contribute by:
    * Reporting bugs or issues
    * Suggesting improvements to the code
    * Contributing new features or analyses

**8. License**

* This code is released under the MIT License

**9. Contact**

* For any inquiries, please contact me at hi at pablohernandezb.dev
