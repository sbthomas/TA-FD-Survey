CC = cc
CFLAGS = -Wall
LDFLAGS = -lm
SRC = survey_processor.c

all: survey_summary survey_distances survey_full survey_averager

survey_summary: $(SRC)
	$(CC) $(CFLAGS) -DSUMMARY -o $@ $< $(LDFLAGS)

survey_distances: $(SRC)
	$(CC) $(CFLAGS) -DDISTANCES -o $@ $< $(LDFLAGS)

survey_full: $(SRC)
	$(CC) $(CFLAGS) -DFULL -o $@ $< $(LDFLAGS)

survey_averager: survey_averager.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f survey_summary survey_distances survey_full survey_averager

.PHONY: all clean
