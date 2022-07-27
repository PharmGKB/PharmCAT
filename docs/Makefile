setup:
	bundle config set --local path 'vendor/bundle'
	bundle install

serve:
	bundle exec jekyll serve --incremental

build:
	bundle exec jekyll build

clean:
	rm -rf _site

next: clean build
	surge _site next.pharmcat.org
