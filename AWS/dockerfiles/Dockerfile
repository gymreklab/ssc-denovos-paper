FROM ileena/cookiemonstr:latest

# Add the fetch_and_run.sh script
ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh

# Get set up to run
WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
