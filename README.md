# IGVF CRISPR Jamboree 2024

Welcome to this year's CRISPR Jamboree! We will be using this repository as a central location for submitting your code. You will be provided with a cluster environment for development during the Jamboree. Below are instructions for accessing and setting up your computational environment, as well as submitting your code to this repository.

## Accessing the cluster

1. Access the server at http://34.135.229.230.
2. You will be prompted for your login credentials. For your username, enter the prefix of the email you entered on [this Google Sheet](https://docs.google.com/spreadsheets/d/1u-7joOMWmFn3490ZXBTOmWGfuVvr9sy7jpJ3owXy1Zg/edit#gid=0). For example, for johndoe@gmail.com it is johndoe. You may choose any password.
3. The sample data is located at `/mnt/shared/`.

## Setting up R (if desired)
1. Open a terminal and run `R -e "IRkernel::installspec()"`.
2. Refresh your browser. 

## Establishing GitHub access from the cluster

1. Fork this repository to your personal GitHub account by clicking the `Fork` button [here](https://github.com/IGVF-CRISPR/CRISPR-JAMBOREE).
2. Set up a fine-grained personal access token [here](https://github.com/settings/tokens?type=beta) by clicking `Generate new token`, entering a token name like "IGVF CRISPR Jamboree 2024", setting the expiration to "7 days", clicking `Only select repositories`, selecting the forked `CRISPR-JAMBOREE` repository, clicking `Repository permissions`, scrolling down to `Contents`, clicking `Access: No access`, choosing `Read and write`, and clicking `Generate token`. Leave this page open. You will need to use the token you generated here in step 7 below.
3. Go back to the cluster, and open a terminal.
4. Clone your forked repository via HTTPS, e.g. `git clone https://github.com/ekatsevi/CRISPR-JAMBOREE.git`.
5. Navigate to the repository via `cd CRISPR-JAMBOREE`.
6. Set up `Git` by issuing the following commands:
- `git config --global user.name "Your name"`, replacing "Your name" with your name
- `git config --global user.email "Your email"`, replacing "Your email" with your email
- `git config --global credential.helper 'store --file ~/.my-credentials'`
7. Then, establish access to `CRISPR-JAMBOREE` on GitHub by entering `git push`, followed by your GitHub username, followed by the personal access token you created in step 2. In particular, you will need to copy and paste the personal access token from step 2 into the terminal. Note that for security reasons, the personal access token will not be displayed when you paste it. Just press enter once you have pasted the token.

## Developing and submitting your code 

1. Within the `CRISPR-JAMBOREE` repository, enter the directory corresponding to your group (`inference`, `guide-assignment`, etc).
2. Within that directory, create a further directory with the name of the method you will be working on (e.g. `wilcoxon-test`).
3. Work on your code within this directory during the Jamboree, committing as usual.
4. At the end of the Jamboree, push your work to your work to your forked repository on GitHub via `git push`.
5. Then, submit a pull request by navigating to your forked repository on GitHub, then clicking `Contribute`, then clicking `Open pull request`, then clicking `Create pull request`. 
