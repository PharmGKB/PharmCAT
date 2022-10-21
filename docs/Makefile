setup:
	bundle config set --local path 'vendor/bundle'
	bundle install

serve:
	bundle exec jekyll serve --incremental --config _config.local.yml

build:
	bundle exec just-the-docs rake search:init
	bundle exec jekyll build --config _config.local.yml

clean:
	rm -rf _site

next: clean build
	surge _site next.pharmcat.org
